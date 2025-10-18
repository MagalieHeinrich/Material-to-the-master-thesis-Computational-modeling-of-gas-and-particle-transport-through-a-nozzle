import ezdxf
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import math
temp_dxf_file_name = "C:/Users/Magalie/Downloads/Richtig_rum.dxf"

# Define functions 
def rotate_points_2d_around_z(points, angle_degrees):
    """
    rotates a list of points aound the z axis, aroung (0,0,0).
    angle_degrees: Angle in degrees. Positive values get rotated mathematically positive etc.
    """
    angle_rad = math.radians(angle_degrees)
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)

    rotated_points = []
    for p in points:
        x, y = p[0], p[1]
        # z is not changed
        z = p[2] if len(p) > 2 else 0.0

        new_x = x * cos_a - y * sin_a
        new_y = x * sin_a + y * cos_a
        rotated_points.append([new_x, new_y, z])
    return rotated_points


def process_dxf_with_polylines(dxf_file_path, rotation_angle_degrees=0):
    """
    Reads DXF-file, extracts information about zu LINE- und LWPOLYLINE-entities,
    rotates the points of the polylines and parametrizizes.

    dxf_file_path: Pfad zur DXF-Datei.
    rotation_angle_degrees: Angle in degrees. Positive values get rotated mathematically positive etc.
    """
    try:
        doc = ezdxf.readfile(dxf_file_path)
    except ezdxf.DXFStructureError:
        print(f"Error: {dxf_file_path} file is not DXF or damaged.")
        return
    except IOError:
        print(f"Error: {dxf_file_path} not found.")
        return

    modelspace = doc.modelspace()
    lwpolyline_count = 0

    print(f"Analyze DXF-File: {dxf_file_path}\n")
    if rotation_angle_degrees != 0:
        print(f"Note: All LWPOLYLINE-Points are rotated by {rotation_angle_degrees} degrees.\n")
    with open(dxf_file_path+".csv", "w") as write_file:
        write_file.write("")

    all_polylines_data = {} # Saves files for all polylines

    for entity in modelspace:
        if entity.dxftype() == 'LWPOLYLINE':
            lwpolyline_count += 1
            original_polyline_points = []
            for vertex in entity.get_points():
                original_polyline_points.append([vertex[0], vertex[1], 0.0]) # Z=0, because 2D

            # --- use the function "rotate_points_2d_around_z" ---
            rotated_polyline_points = rotate_points_2d_around_z(original_polyline_points, rotation_angle_degrees)


            print(f"=== LWPOLYLINE {lwpolyline_count} (Handle: {entity.dxf.handle}) ===")
            print(f"Number of points (after rotation): {len(rotated_polyline_points)}")
            print("Points (X, Y, Z)")
            with open(dxf_file_path+".csv", "a") as write_file:
                write_file.write("\n")
                for i, pt in enumerate(rotated_polyline_points):
                   
                    write_file.write(f"{pt[0]:.4f}; {pt[1]:.4f}\n")
            

            # --- Parametrization of the LWPOLYLINE ---
            points_array = np.array(rotated_polyline_points)

            # generate parameter t (summed arc length)
            diffs = np.diff(points_array, axis=0)
            segment_lengths = np.sqrt(np.sum(diffs**2, axis=1))
            t_values = np.insert(np.cumsum(segment_lengths), 0, 0.0)

            if t_values[-1] > 0:
                t_values_normalized = t_values / t_values[-1]
            else:
                t_values_normalized = t_values

        
            # Cubic spline interpolation for X, Y und Z als functions of t
            try:
                cs_x = CubicSpline(t_values_normalized, points_array[:, 0])
                cs_y = CubicSpline(t_values_normalized, points_array[:, 1])
                cs_z = CubicSpline(t_values_normalized, points_array[:, 2])

                # Save the spline functions for later
                all_polylines_data[entity.dxf.handle] = {
                    'original_points': original_polyline_points,
                    'rotated_points': points_array,
                    't_values': t_values_normalized,
                    'cs_x': cs_x,
                    'cs_y': cs_y,
                    'cs_z': cs_z
                }



            except Exception as e:
                print(f"  Error at spline interpolation for the LWPOLYLINE: {e}\n")

    if lwpolyline_count == 0:
        print("Did not find any LWPOLYLINE-entities in the DXF-file.")

    return all_polylines_data

# Main Part
param_data = process_dxf_with_polylines(temp_dxf_file_name, rotation_angle_degrees=-90)

# Optional visualization of the parameterized curves (needs matplotlib) ---
if param_data:
    plt.figure(figsize=(10, 8))
    for handle, data in param_data.items():
        # Plotting of the rotated points and the interpolation of the rotated points
        rotated_points_array = data['rotated_points']
        cs_x = data['cs_x']
        cs_y = data['cs_y']
        
        plt.plot(rotated_points_array[:, 0], rotated_points_array[:, 1], 'o', label=f'LWPolyline {handle} (rotated points)')

        # Plotting of the evaluated points of the spline interpolation
        t_fine = np.linspace(0, 1, 500)
        x_fine = cs_x(t_fine)
        y_fine = cs_y(t_fine)
        plt.plot(x_fine, y_fine, '-', label=f'LWPolyline {handle} (Spline-Interpolation)')

    plt.xlabel('X-coordinate (after rotation)')
    plt.ylabel('Y-coordinate (after rotation)')
    plt.title('Rotated LWPOLYLINEs from DXF and their spline interpolation')
    plt.grid(True)
    plt.legend()
    plt.axis('equal') # Important for plotting the correct aspect ratio
    plt.show()