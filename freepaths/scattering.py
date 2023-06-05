"""
Modules provides scattering processes on various complex objects.
Each function determines whether the scattering should happen and
call corresponding function for scattering on corresponding primitive.
"""

from math import pi, cos, sin, tan, exp, sqrt, atan, asin, acos
from random import random
from numpy import sign

from freepaths.config import cf
from freepaths.move import move
from freepaths.scattering_types import Scattering
from freepaths.scattering_primitives import *
from freepaths.scattering_parabolic import *
from freepaths.scatterers import *


def internal_scattering(ph, flight, scattering_types):
    """Check if the time passed since previous diffuse scattering event reached
    the time until an internal scattering event, and if yes, scatters randomly"""
    if flight.time_since_previous_scattering >= ph.time_of_internal_scattering:
        scattering_types.internal = random_scattering(ph)


def reinitialization(ph, scattering_types):
    """Re-thermalize (diffusely) phonon when it comes back to one of the hot sides"""
    x, y, _ = move(ph, cf.timestep)

    if cf.hot_side_position_bottom and y < 0:
        scattering_types.hot_side = horizontal_surface_up_scattering(ph, cf.side_wall_roughness, is_diffuse=True)

    if cf.hot_side_position_top and y > cf.length:
        scattering_types.hot_side = horizontal_surface_down_scattering(ph, cf.side_wall_roughness, is_diffuse=True)

    if cf.hot_side_position_right and x > cf.width/2:
        scattering_types.hot_side = vertical_surface_left_scattering(ph, cf.side_wall_roughness, is_diffuse=True)

    if cf.hot_side_position_left and x < - cf.width/2:
        scattering_types.hot_side = vertical_surface_right_scattering(ph, cf.side_wall_roughness, is_diffuse=True)


def scattering_on_circular_holes(ph, x0, y0, radius, scattering_types, x, y, z):
    """Check if a phonon strikes a circular hole and calculate the new direction"""

    # If phonon is inside the circle with given radius:
    if (x - x0)**2 + (y - y0)**2 <= radius**2:
        if y == y0: y += 1e-9 # Prevent division by zero
        tangent_theta = atan((x - x0)/(y - y0))
        scattering_types.holes = circle_outer_scattering(ph, tangent_theta, y, y0, cf.hole_roughness)

def scattering_on_semicircular_holes(ph, x0, y0, R, scattering_types, x, y, z):
    """Check if a phonon strikes a circular hole and calculate the new direction"""

    # If phonon is inside the circle with radius R:
    if (x - x0)**2 + (y - y0)**2 <= R**2 and x>=x0:
        r1=abs(y0 - y + (x-x0)*cos(ph.theta)/abs(sin(ph.theta)))
        if R>r1:
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*sin(abs(ph.theta)))  # Angle to the normal to the surface
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                scattering_types.holes = Scattering.SPECULAR
                ph.theta = - ph.theta

            # Diffuse scattering:
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1
                    # Lambert cosine distribution:
                    ph.theta = - sign(sin(ph.theta))*pi/2 + asin(2*random()-1)
                    ph.phi = asin((asin(2*random() - 1))/(pi/2))

                    # Accept the angles only if they do not lead to new scattering:
                    if no_new_scattering(ph):
                        break
        else:
            # Calculate angle to the surface and specular scattering probability:
            if y == y0: y += 1e-9 # Prevent division by zero
            tangent_theta = atan((x - x0)/(y - y0))
            a = acos(cos(ph.phi)*cos(ph.theta + sign(y - y0)*tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:
                ph.theta = - ph.theta - pi + 2*tangent_theta
                scattering_types.holes = Scattering.SPECULAR

            # Diffuse scattering:
            else:
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1

                    # Random distribution:
                    # theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    # phi = asin(2*random() - 1)

                    # Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
                    ph.phi = asin((asin(2*random() - 1))/(pi/2))

                    # Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph):
                        break
def scattering_on_arccircular_v_holes(ph, x0, y0, R ,Rinner,alphap, scattering_types, x, y, z):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:# to prevent division by 0 
        x = x + 1e-12
    theta0 = atan((y-y0)/(x-x0))
    if y == y0:# to prevent division by 0 
        y = y + 1e-12
    tangent_theta = atan((x-x0)/(y-y0)) # use for scatering sepcular of circular boundary
    
    #to know if it is in or not
    if  (Rinner**2 <=(x - x0)**2 + (y - y0)**2 <= R**2) and x>=x0 and -alphap/2 <= theta0 <= alphap/2:
        xp=ph.x-x0
        yp=ph.y-y0
        thetapre = atan((yp)/(xp)) #angle of previous point 
 
        
        #Parameter
        continu = 0 # to know if already scatered
        beta= pi/2-alphap/2 #angle for triangle
        pillar_wall_angle = pi/2 
        
       
        
        if  yp > 0  and (R**2 <= xp**2 + yp**2 and -alphap/2 <= thetapre <= alphap/2)==False and (Rinner**2 >= xp**2 + yp**2 and -alphap/2  < thetapre <  alphap/2)==False:
            continu=1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta -(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
            
            # Specular scattering:
            if random() <p:
                
                ph.theta = (- ph.theta +2*beta)
                scattering_types.holes = Scattering.SPECULAR
                
            else: 
                ph.theta = asin(2*random() - 1) -1*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = Scattering.DIFFUSE
           
                
       
        if  yp < 0  and (R**2 <= xp**2 + yp**2 and -alphap/2 <= thetapre <= alphap/2)==False and(Rinner**2 >= xp**2 + yp**2 and -alphap/2 < thetapre <  alphap/2)==False:
            continu=1
            
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta +(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random() <p:
                ph.theta = - ph.theta - 2*beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                #rand_sign = sign((2*random() - 1))
                ph.theta = pi - asin(random()) +1*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes= Scattering.DIFFUSE

            
        
        if R**2 <= xp**2 + yp**2 and continu == 0: 
           
            
            continu=1
            if y == y0: 
                y += 1e-9 # Prevent division by zero
            tangent_theta = atan((x - x0)/(y - y0))
            a = acos(cos(ph.phi)*cos(ph.theta + sign(y - y0)*tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random() <p:
                ph.theta = - ph.theta - pi + 2*tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else: 
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1
        
                     #Random distribution:
                    #theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    #phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
                    ph.phi = asin((asin(2*random() - 1))/(pi/2))
                        
                     #Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph):
                         break
      
        if continu ==0:
            tangent_theta = atan((x - x0)/(y - y0))
            a = atan(tan((pi/2 - ph.theta) + tangent_theta) * cos(ph.phi - (pi / 2 - pillar_wall_angle)))
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:

                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0))**2 + (abs(y) - abs(y0))**2) >= sqrt((abs(ph.x) - abs(x0))**2 + (abs(ph.y) - abs(y0))**2) :

                    # If theta does not reflect back:
                    #if ph.phi < pi/2 - 2 * pillar_wall_angle:
                        #ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    #else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR
            
            else: 
                #attempt = 0
                #while attempt < 10:
                    #attempt += 1
                   
               
                ph.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
                ph.phi = asin((asin(2*random() - 1))/(pi/2)) - (pi / 2 - cf.pillar_wall_angle)
                scattering_types.pillars = Scattering.DIFFUSE
                    #if no_new_scattering(ph):
                         #break

def scattering_on_arccircular_v_demi_down_holes(ph, x0, y0, R ,Rinner,alphap,alphap2, scattering_types, x, y, z):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:# to prevent division by 0 
        x = x + 1e-12
    theta0 = atan((y-y0)/(x-x0))
    if y == y0:# to prevent division by 0 
        y = y + 1e-12
    tangent_theta = atan((x-x0)/(y-y0)) # use for scatering sepcular of circular boundary
    
    #to know if it is in or not
    if  (Rinner**2 <=(x - x0)**2 + (y - y0)**2 <= R**2) and x>=x0 and -alphap/2 <= theta0 <= -alphap2/2:
        xp=ph.x-x0
        yp=ph.y-y0
        thetapre = atan((yp)/(xp)) #angle of previous point 
      
        #Parameter
        continu = 0 # to know if already scatered
        beta= pi/2-alphap/2 #angle for triangle
        beta2=pi/2- alphap2/2
        pillar_wall_angle = pi/2 
        y_mid= -xp*tan(alphap/2)+10e-9#-alphap2/2)
       
        
        if  yp > y_mid  and (R**2 <= xp**2 + yp**2 and -alphap/2 <= thetapre <= -alphap2/2)==False and (Rinner**2 >= xp**2 + yp**2 and -alphap/2  < thetapre <  -alphap2/2)==False:
            continu=1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta +(pi/2 - beta2)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
            
            # Specular scattering:
            if random()<p:
                
                ph.theta =(- ph.theta -2*beta2)
                scattering_types.holes = Scattering.SPECULAR
                
            else: 
                ph.theta = asin(2*random() - 1) +1*(pi/2 - beta2)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = Scattering.DIFFUSE
           
                
       
        if  yp < y_mid  and (R**2 <= xp**2 + yp**2 and -alphap/2 <= thetapre <= -alphap2/2)==False and(Rinner**2 >= xp**2 + yp**2 and -alphap/2 < thetapre <  -alphap2/2)==False:
            continu=1
            
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta +(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random() <p:
                ph.theta = - ph.theta - 2*beta
                scattering_types.holes = Scattering.SPECULAR
            else:
                #rand_sign = sign((2*random() - 1))
                ph.theta = pi - asin(random()) +1*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes= Scattering.DIFFUSE

            
        
        if R**2 <= xp**2 + yp**2 and continu == 0: 
            #for sol in S:
                #Ima_sol = sol- conj(sol)
                #if Ima_sol==0:
                    #if -1<=sol<=1:
                        #if abs(acos(sol))<=alphap/2:
            
            continu=1
            if y == y0: 
                y += 1e-9 # Prevent division by zero
            tangent_theta = atan((x - x0)/(y - y0))
            a = acos(cos(ph.phi)*cos(ph.theta + sign(y - y0)*tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random() <p:
                ph.theta = - ph.theta - pi + 2*tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else: 
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1
        
                     #Random distribution:
                    #theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    #phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
                    ph.phi = asin((asin(2*random() - 1))/(pi/2))
                        
                     #Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph):
                         break
      
        if continu ==0:
            tangent_theta = atan((x - x0)/(y - y0))
            a = atan(tan((pi/2 - ph.theta) + tangent_theta) * cos(ph.phi - (pi / 2 - pillar_wall_angle)))
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random() <p:

                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0))**2 + (abs(y) - abs(y0))**2) >= sqrt((abs(ph.x) - abs(x0))**2 + (abs(ph.y) - abs(y0))**2) :

                    # If theta does not reflect back:
                    #if ph.phi < pi/2 - 2 * pillar_wall_angle:
                        #ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    #else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR
            
            else: 
                #attempt = 0
                #while attempt < 10:
                    #attempt += 1
                   
               
                ph.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
                ph.phi = asin((asin(2*random() - 1))/(pi/2)) - (pi / 2 - cf.pillar_wall_angle)
                scattering_types.pillars = Scattering.DIFFUSE
                    #if no_new_scattering(ph):
                        
def scattering_on_arccircular_v_demi_up_holes(ph, x0, y0, R ,Rinner,alphap,alphap2, scattering_types, x, y, z):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:# to prevent division by 0 
        x = x + 1e-12
    theta0 = atan((y-y0)/(x-x0))
    if y == y0:# to prevent division by 0 
        y = y + 1e-12
    tangent_theta = atan((x-x0)/(y-y0)) # use for scatering sepcular of circular boundary
    
    #to know if it is in or not
    if  (Rinner**2 <=(x - x0)**2 + (y - y0)**2 <= R**2) and x>=x0 and alphap2/2 <= theta0 <= alphap/2:
        xp=ph.x-x0
        yp=ph.y-y0
        thetapre = atan((yp)/(xp)) #angle of previous point 
      
        #Parameter
        continu = 0 # to know if already scatered
        beta= pi/2-alphap/2 #angle for triangle
        beta2= pi/2-alphap2/2
        pillar_wall_angle = pi/2 
        y_mid= xp*tan(alphap/2)-10e-9#-alphap2/2)
       
        
        if  yp > y_mid  and (R**2 <= xp**2 + yp**2 and alphap2/2 <= thetapre <= alphap/2)==False and (Rinner**2 >= xp**2 + yp**2 and alphap2/2  < thetapre <  alphap/2)==False:
            continu=1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta -(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
            
            # Specular scattering:
            if random() <p:
                
                ph.theta =(- ph.theta +2*beta)
                scattering_types.holes = Scattering.SPECULAR
                
            else: 
                ph.theta = asin(2*random() - 1) -1*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = Scattering.DIFFUSE
           
                
       
        if  yp < y_mid  and (R**2 <= xp**2 + yp**2 and alphap2/2 <= thetapre <= alphap/2)==False and(Rinner**2 >= xp**2 + yp**2 and alphap2/2 < thetapre <  alphap/2)==False:
            continu=1
            
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta -(pi/2 - beta2)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random() <p:
                ph.theta = - ph.theta + 2*beta2
                scattering_types.holes = Scattering.SPECULAR
            else:
                #rand_sign = sign((2*random() - 1))
                ph.theta = pi - asin(random()) -1*(pi/2 - beta2)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes= Scattering.DIFFUSE

            
        
        if R**2 <= xp**2 + yp**2 and continu == 0: 
            #for sol in S:
                #Ima_sol = sol- conj(sol)
                #if Ima_sol==0:
                    #if -1<=sol<=1:
                        #if abs(acos(sol))<=alphap/2:
            
            continu=1
            if y == y0: 
                y += 1e-9 # Prevent division by zero
            tangent_theta = atan((x - x0)/(y - y0))
            a = acos(cos(ph.phi)*cos(ph.theta + sign(y - y0)*tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random() <p:
                ph.theta = - ph.theta - pi + 2*tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else: 
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1
        
                     #Random distribution:
                    #theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    #phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
                    ph.phi = asin((asin(2*random() - 1))/(pi/2))
                        
                     #Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph):
                         break
      
        if continu ==0:
            tangent_theta = atan((x - x0)/(y - y0))
            a = atan(tan((pi/2 - ph.theta) + tangent_theta) * cos(ph.phi - (pi / 2 - pillar_wall_angle)))
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random()< p:

                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0))**2 + (abs(y) - abs(y0))**2) >= sqrt((abs(ph.x) - abs(x0))**2 + (abs(ph.y) - abs(y0))**2) :

                    # If theta does not reflect back:
                    #if ph.phi < pi/2 - 2 * pillar_wall_angle:
                        #ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    #else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR
            
            else: 
                #attempt = 0
                #while attempt < 10:
                    #attempt += 1
                   
               
                ph.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
                ph.phi = asin((asin(2*random() - 1))/(pi/2)) - (pi / 2 - cf.pillar_wall_angle)
                scattering_types.pillars = Scattering.DIFFUSE
                    #if no_new_scattering(ph):
                         #break                                 #break        
def scattering_on_arccircular_h_holes(ph, x0, y0, R ,Rinner,alphap, scattering_types, x, y, z):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:# to prevent division by 0 
        x = x + 1e-12
    theta0 = atan((x-x0)/(y-y0))
    if y == y0:# to prevent division by 0 
        y = y + 1e-12
    tangent_theta = atan((x-x0)/(y-y0))
    
    #to know if it is in or not
    if  (Rinner**2 <=(x - x0)**2 + (y - y0)**2 <= R**2) and y>=y0 and -alphap/2 <= theta0 <= alphap/2:
        xp=ph.x-x0
        yp=ph.y-y0
        thetapre = atan((xp)/(yp))
       
    
        continu = 0
        beta= alphap/2
        pillar_wall_angle = pi/2
        
       
        
        if  xp > 0  and (R**2 <= xp**2 + yp**2 and -alphap/2 <= thetapre <= alphap/2)==False and (Rinner**2 >= xp**2 + yp**2 and  -alphap/2  < thetapre <  alphap/2)==False:
            continu=1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta -(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
            
            # Specular scattering:
            if random()<p:
                
                ph.theta = (- ph.theta +2*beta)
                scattering_types.holes = Scattering.SPECULAR
                
            else: 
                ph.theta = pi-(asin(2*random() - 1) +1*(pi/2 - beta))
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = Scattering.DIFFUSE
           
                
       
        if  xp < 0  and (R**2 <= xp**2 + yp**2 and -alphap/2 <= thetapre <= alphap/2)==False and(Rinner**2 >= xp**2 + yp**2 and -alphap/2 < thetapre <  alphap/2)==False:
            continu=1
            
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta +(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random()<p:
                ph.theta = (-ph.theta - 2*beta)
                scattering_types.holes = Scattering.SPECULAR
            else:
                #rand_sign = sign((2*random() - 1))
                ph.theta = pi - asin(random()) +1*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes= Scattering.DIFFUSE

            
        
        if R**2 <= xp**2 + yp**2 and continu == 0: 
            #for sol in S:
                #Ima_sol = sol- conj(sol)
                #if Ima_sol==0:
                    #if -1<=sol<=1:
                        #if abs(acos(sol))<=alphap/2:
            
            continu=1
            if y == y0: 
                y += 1e-9 # Prevent division by zero
            tangent_theta = atan((x - x0)/(y - y0))
            a = acos(cos(ph.phi)*cos(ph.theta + sign(y - y0)*tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random() <p:
                ph.theta = - ph.theta - pi + 2*tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else: 
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1
        
                     #Random distribution:
                    #theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    #phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
                    ph.phi = asin((asin(2*random() - 1))/(pi/2))
                        
                     #Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph):
                         break
      
        if continu ==0:
            tangent_theta = atan((x - x0)/(y - y0))
            a = atan(tan((pi/2 - ph.theta) + tangent_theta) * cos(ph.phi - (pi / 2 - pillar_wall_angle)))
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random() < p:

                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0))**2 + (abs(y) - abs(y0))**2) >= sqrt((abs(ph.x) - abs(x0))**2 + (abs(ph.y) - abs(y0))**2) :

                    # If theta does not reflect back:
                    #if ph.phi < pi/2 - 2 * pillar_wall_angle:
                        #ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    #else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR
            
            else: 
                #attempt = 0
                #while attempt < 10:
                    #attempt += 1
                   
               
                ph.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
                ph.phi = asin((asin(2*random() - 1))/(pi/2)) - (pi / 2 - cf.pillar_wall_angle)
                scattering_types.pillars = Scattering.DIFFUSE
                    #if no_new_scattering(ph):
                        
def scattering_on_arccircular_h_reverse_holes(ph, x0, y0, R ,Rinner,alphap, scattering_types, x, y, z):
    """Check if a phonon strikes a circular hole and calculate the new direction"""
    if x == x0:# to prevent division by 0 
        x = x + 1e-12
    theta0 = atan((x-x0)/(y-y0))
    if y == y0:# to prevent division by 0 
        y = y + 1e-12
    tangent_theta = atan((x-x0)/(y-y0))
    
    #to know if it is in or not
    if  (Rinner**2 <=(x - x0)**2 + (y - y0)**2 <= R**2) and y<=y0 and -alphap/2 <= theta0 <= alphap/2:
        xp=ph.x-x0
        yp=ph.y-y0
        thetapre = atan((xp)/(yp))
       
    
        continu = 0
        beta= alphap/2
        pillar_wall_angle = pi/2
        
       
        
        if  xp > 0  and (R**2 <= xp**2 + yp**2 and -alphap/2 <= thetapre <= alphap/2)==False and (Rinner**2 >= xp**2 + yp**2 and  -alphap/2  < thetapre <  alphap/2)==False:
            continu=1
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta +(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
            
            # Specular scattering:
            if random()<p:
                
                ph.theta = (- ph.theta -2*beta)
                scattering_types.holes = Scattering.SPECULAR
                
            else: 
                ph.theta = asin(2*random() - 1) +1*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = Scattering.DIFFUSE
           
                
       
        if  xp < 0  and (R**2 <= xp**2 + yp**2 and -alphap/2 <= thetapre <= alphap/2)==False and(Rinner**2 >= xp**2 + yp**2 and -alphap/2 < thetapre <  alphap/2)==False:
            continu=1
            
            # Calculate angle to the surface and specular scattering probability:
            a = acos(cos(ph.phi)*cos(ph.theta -(pi/2 - beta)))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random()<p:
                ph.theta = (-ph.theta + 2*beta)
                scattering_types.holes = Scattering.SPECULAR
            else:
                #rand_sign = sign((2*random() - 1))
                
                ph.theta = asin(2*random() - 1) -1*(pi/2 - beta)
                ph.phi = asin((asin(2*random() - 1))/(pi/2))
                scattering_types.holes = Scattering.DIFFUSE

            
        
        if R**2 <= xp**2 + yp**2 and continu == 0: 
            #for sol in S:
                #Ima_sol = sol- conj(sol)
                #if Ima_sol==0:
                    #if -1<=sol<=1:
                        #if abs(acos(sol))<=alphap/2:
            
            continu=1
            if y == y0: 
                y += 1e-9 # Prevent division by zero
            tangent_theta = atan((x - x0)/(y - y0))
            a = acos(cos(ph.phi)*cos(ph.theta + sign(y - y0)*tangent_theta))
            p = specularity(a, cf.hole_roughness, ph.wavelength)
 
            # Specular scattering:
            if random() <p:
                ph.theta = - ph.theta - pi + 2*tangent_theta
                scattering_types.holes = Scattering.SPECULAR
            else: 
                scattering_types.holes = Scattering.DIFFUSE
                attempt = 0
                while attempt < 10:
                    attempt += 1
        
                     #Random distribution:
                    #theta = tangent_theta - (-pi/2 + pi*random()) - pi*(y < y0)
                    #phi = asin(2*random() - 1)
                    # #Lambert cosine distribution:
                    ph.theta = tangent_theta - (asin(2*random() - 1)) - pi*(y < y0)
                    ph.phi = asin((asin(2*random() - 1))/(pi/2))
                        
                     #Accept the angles only if they do not immediately cause new scattering:
                    if no_new_scattering(ph):
                         break
      
        if continu ==0:
            tangent_theta = atan((x - x0)/(y - y0))
            a = atan(tan((pi/2 - ph.theta) + tangent_theta) * cos(ph.phi - (pi / 2 - pillar_wall_angle)))
            p = specularity(a, cf.pillar_roughness, ph.wavelength)

            # Specular scattering:
            if random()<p:

                # If phonon moves from the center of the pillar to the wall:
                if sqrt((abs(x) - abs(x0))**2 + (abs(y) - abs(y0))**2) >= sqrt((abs(ph.x) - abs(x0))**2 + (abs(ph.y) - abs(y0))**2) :

                    # If theta does not reflect back:
                    #if ph.phi < pi/2 - 2 * pillar_wall_angle:
                        #ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                    # Regular reflection:
                    #else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = ph.phi - (pi / 2 - pillar_wall_angle)

                # If phonon strikes the wall as it goes towards the center:
                else:
                    ph.theta = - ph.theta - pi + 2*tangent_theta
                    ph.phi = -sign(ph.phi) * ph.phi - 2 * pillar_wall_angle
                scattering_types.holes = Scattering.SPECULAR
            
            else: 
                #attempt = 0
                #while attempt < 10:
                    #attempt += 1
                   
               
                ph.theta = tangent_theta - asin(2*random()-1) + pi*(y >= y0)
                ph.phi = asin((asin(2*random() - 1))/(pi/2)) - (pi / 2 - cf.pillar_wall_angle)
                scattering_types.pillars = Scattering.DIFFUSE
                    #if no_new_scattering(ph):
                         #break 
                         
def scattering_on_rectangular_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z):
    """Check if the phonon strikes a rectangular hole and calculate new direction"""

    # If the phonon is inside the rectangle:
    if (abs(x - x0) <= Lx / 2) and (abs(y - y0) <= Ly / 2):

        # Coordinate y of the intersection with the hole side:
        y1 = (y0 - y) + cos(ph.theta)*(Lx/2 - abs(x0 - x))/abs(sin(ph.theta))

        # Scattering on the left wall:
        if abs(y1) <= Ly/2 and x < x0:
            scattering_types.holes = vertical_surface_left_scattering(ph, cf.hole_roughness)

        # Scattering on the right wall:
        elif abs(y1) <= Ly/2 and x > x0:
            scattering_types.holes = vertical_surface_right_scattering(ph, cf.hole_roughness)

        # Scattering on the top wall:
        elif y > y0:
            scattering_types.holes = horizontal_surface_up_scattering(ph, cf.side_wall_roughness)

        # Scattering on the bottom wall:
        else:
            scattering_types.holes = horizontal_surface_down_scattering(ph, cf.side_wall_roughness)


def scattering_on_circular_pillars(ph, pillar, scattering_types, x, y, z):
    """Check if a phonon strikes a circular pillar and calculate new direction"""

    # Cone radius at a given z coordinate:
    radius = pillar.diameter/2# - (z - cf.thickness / 2) / tan(pillar.wall_angle)
    distance_from_pillar_center = sqrt((x - pillar.x)**2 + (y - pillar.y)**2)
    distance_from_pillar_center_original = sqrt((ph.x - pillar.x)**2 + (ph.y - pillar.y)**2)
    step = 2 * ph.speed * cf.timestep

    # If phonon crosses the pillar boundary. Third condition is to exclude all other pillars:
    if (distance_from_pillar_center >= radius and z > cf.thickness / 2 and distance_from_pillar_center < radius + step):

        # Calculate angle to the surface and specular scattering probability:
        tangent_theta = atan((x - pillar.x)/(y - pillar.y))
        scattering_types.pillars = circle_inner_scattering(ph, tangent_theta, y, pillar.y, cf.pillar_roughness)


def scattering_on_triangle_down_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction after the scattering"""

    # Angle of the triangle:
    beta = atan(0.5 * Lx / Ly)

    # If phonon is inside the triangle:
    if (Ly/2 - (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):

        # Scattering on the top wall of the triangle:
        if (ph.y > y0 + Ly / 2) and (abs(ph.theta) > pi / 2):
            scattering_types.holes = horizontal_surface_up_scattering(ph, cf.hole_roughness)

        # Scattering on the sidewalls of the triangle:
        else:
            scattering_types.holes = inclined_surfaces_down_scattering(ph, beta, x, x0, cf.hole_roughness)


def scattering_on_triangle_up_holes(ph, x0, y0, Lx, Ly, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

    # Angle of the triangle:
    beta = atan(0.5 * Lx / Ly)

    # If phonon is inside the triangle:
    if (Ly/2 + (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2):

        # Scattering on the bottom wall of the triangle:
        if (ph.y < y0 - Ly / 2) and (abs(ph.theta) < pi / 2):
            scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)

        # Scattering on the sidewalls of the triangle:
        else:
            scattering_types.holes = inclined_surfaces_up_scattering(ph, beta, x, x0, cf.hole_roughness)


def scattering_on_half_triangle_up_holes(ph, x0, y0, Lx, Ly, is_right_half, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

    # Angle of the triangle:
    beta = atan(0.5 * Lx / Ly)

    # If phonon is inside the right side of the triangle:
    if is_right_half and (Ly/2 + (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2) and (x > x0):

            # Scattering on the bottom wall of the triangle:
            if (ph.y < y0 - Ly / 2) and (abs(ph.theta) < pi / 2):
                scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x < x0:
                scattering_types.holes = vertical_surface_left_scattering(ph, cf.hole_roughness)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_up_scattering(ph, beta, x, x0, cf.hole_roughness)

    # If phonon is inside the left side of the triangle:
    if not is_right_half and (Ly/2 + (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2) and (x < x0):

            # Scattering on the bottom wall of the triangle:
            if (ph.y < y0 - Ly / 2) and (abs(ph.theta) < pi / 2):
                scattering_types.holes = horizontal_surface_down_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x > x0:
                scattering_types.holes = vertical_surface_right_scattering(ph, cf.hole_roughness)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_up_scattering(ph, beta, x, x0, cf.hole_roughness)


def scattering_on_half_triangle_down_holes(ph, x0, y0, Lx, Ly, is_right_half, scattering_types, x, y, z):
    """Check if the phonon strikes a reverse triangular hole and calculate new direction"""

    # Angle of the triangle:
    beta = atan(0.5 * Lx / Ly)

    # If phonon is inside the right side of the triangle:
    if is_right_half and (Ly/2 - (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2) and (x > x0):

            # Scattering on the top wall of the triangle:
            if (ph.y > y0 + Ly / 2) and (abs(ph.theta) > pi / 2):
                scattering_types.holes = horizontal_surface_up_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x < x0:
                scattering_types.holes = vertical_surface_left_scattering(ph, cf.hole_roughness)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_down_scattering(ph, beta, x, x0, cf.hole_roughness)

    # If phonon is inside the left side of the triangle:
    if not is_right_half and (Ly/2 - (y - y0) <= (Lx/2 - abs(x - x0))/tan(beta)) and (abs(y - y0) < Ly/2) and (x < x0):

            # Scattering on the top wall of the triangle:
            if (ph.y > y0 + Ly / 2) and (abs(ph.theta) > pi / 2):
                scattering_types.holes = horizontal_surface_up_scattering(ph, cf.hole_roughness)

            # Scattering on the vertical sidewall of the triangle:
            elif ph.x > x0:
                scattering_types.holes = vertical_surface_right_scattering(ph, cf.hole_roughness)

            # Scattering on the inclined sidewall of the triangle:
            else:
                scattering_types.holes = inclined_surfaces_down_scattering(ph, beta, x, x0, cf.hole_roughness)


def scattering_on_right_sidewall(ph, scattering_types, x, y, z):
    """Scatter phonon if it reached right side wall"""
    if x > cf.width/2:
        scattering_types.walls = vertical_surface_left_scattering(ph, cf.side_wall_roughness)


def scattering_on_left_sidewall(ph, scattering_types, x, y, z):
    """Scatter phonon if it reached left side wall"""
    if x < -cf.width/2:
        scattering_types.walls = vertical_surface_right_scattering(ph, cf.side_wall_roughness)


def scattering_on_top_sidewall(ph, scattering_types, x, y, z):
    """Check if the phonon hits top side wall and output new vector"""
    if y > cf.length:
        scattering_types.walls = horizontal_surface_down_scattering(ph, cf.side_wall_roughness)


def scattering_on_bottom_sidewall(ph, scattering_types, x, y, z):
    """Check if the phonon hits bottom side wall and output new vector"""
    if y < 0.0:
        scattering_types.walls = horizontal_surface_up_scattering(ph, cf.side_wall_roughness)


def floor_scattering(ph, scattering_types, x, y, z):
    """Check if the phonon hits the floor surface and calculate new angles"""
    if z < -cf.thickness/2:
        scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)


def ceiling_scattering(ph, scattering_types, x, y, z):
    """Check if the phonon hits the ceiling surface and if this place has a pillar and output new vector"""
    if z > cf.thickness / 2:

        # Regular scattering if there are no pillars:
        if not cf.pillars:
            scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)

        # If it is under a pillar, then scatter from the pillar's ceiling:
        else:
            for pillar in cf.pillars:
                distance_from_pillar_center = sqrt((x - pillar.x)**2 + (y - pillar.y)**2)
                is_under_pillar = (distance_from_pillar_center < pillar.diameter/2)
                if is_under_pillar and z > pillar.height + cf.thickness / 2:
                    scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)
                    return
                elif is_under_pillar and z <= pillar.height + cf.thickness / 2:
                    return
                else:
                    pass

            # Regular scattering if phonon is not under any of the pillars:
            if ph.z < cf.thickness/2:
                scattering_types.top_bottom = in_plane_surface_scattering(ph, cf.top_roughness)


def surface_scattering(ph, scattering_types):
    """Check for a surface scattering on this step"""

    # Preliminary move to see if phonon would cross something:
    x, y, z = move(ph, cf.timestep)

    # Scattering on top and bottom surfaces:
    ceiling_scattering(ph, scattering_types, x, y, z)
    floor_scattering(ph, scattering_types, x, y, z)

    # Scattering on sidewalls:
    if cf.include_right_sidewall:
        scattering_on_right_sidewall(ph, scattering_types, x, y, z)
    if cf.include_left_sidewall:
        scattering_on_left_sidewall(ph, scattering_types, x, y, z)
    if cf.include_top_sidewall:
        scattering_on_top_sidewall(ph, scattering_types, x, y, z)
    if cf.include_bottom_sidewall:
        scattering_on_bottom_sidewall(ph, scattering_types, x, y, z)


    # Scattering on holes:
    if cf.holes:
        # Check for each hole and each hole type:
        for hole in cf.holes:
            if isinstance(hole, CircularHole):
                scattering_on_circular_holes(ph, hole.x, hole.y, hole.diameter/2, scattering_types, x, y, z)

            elif isinstance(hole, RectangularHole):
                scattering_on_rectangular_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, scattering_types, x, y, z)

            elif isinstance(hole, TriangularUpHole):
                scattering_on_triangle_up_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, scattering_types, x, y, z)

            elif isinstance(hole, TriangularDownHole):
                scattering_on_triangle_down_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, scattering_types, x, y, z)

            elif isinstance(hole, TriangularUpHalfHole):
                scattering_on_half_triangle_up_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, hole.is_right_half, scattering_types, x, y, z)

            elif isinstance(hole, TriangularDownHalfHole):
                scattering_on_half_triangle_down_holes(ph, hole.x, hole.y, hole.size_x, hole.size_y, hole.is_right_half, scattering_types, x, y, z)

            elif isinstance(hole, ParabolaTop):
                top_parabola_scattering(ph, hole, cf.side_wall_roughness, scattering_types, x, y, z)

            elif isinstance(hole, ParabolaBottom):
                bottom_parabola_scattering(ph, hole, cf.side_wall_roughness, scattering_types, x, y, z)

            else:
                pass

            # If there was any scattering, then no need to check rest of the holes:
            if scattering_types.holes is not None:
                break

    # Check for each pillar:
    if cf.pillars:
        for pillar in cf.pillars:
            if isinstance(pillar, CircularPillar):
                scattering_on_circular_pillars(ph, pillar, scattering_types, x, y, z)

            # If there was any scattering, then no need to check other pillars:
            if scattering_types.pillars is not None:
                break

    # Correct angle if it became more than 180 degrees:
    ph.correct_angle()
