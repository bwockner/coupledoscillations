import csv
import bpy
from math import sin, cos, pi

scene   = bpy.context.scene

track = bpy.data.objects['Circle']
ball = bpy.data.objects['Sphere']

with open('D:\\Uni Work\\MECH3200\\DATA.txt', 'r') as f:
    reader = csv.reader(f)
    reader2 = list(reader)
        
time = [float(x) for x in reader2[0]]
x_disp = [float(x) for x in reader2[1]]
theta = [float(x) for x in reader2[2]]

R_t = 1
R_bg = 0.8

scene.frame_start = 1
scene.frame_end = len(time)

s = bpy.data.scenes[0] 

for i in range(len(time)):
    
    bpy.ops.object.select_all(action='DESELECT')
    track.select = True
    ball.select = True
                  
    track.location = (0,-x_disp[i],1)
    track.rotation_euler = (-pi/2, pi/2 + x_disp[i]/R_t, -pi/2)
    trackloc = track.location
    ball.location = (trackloc.x,trackloc.y - R_bg*sin(theta[i]),R_t-R_bg*cos(theta[i])))
    
    track.keyframe_insert(data_path="location", frame=i)
    track.keyframe_insert(data_path="rotation_euler", frame=i)
    ball.keyframe_insert(data_path="location", frame=i)
    
    
