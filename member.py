''' 
Todo:
 - Welds
 - Plastic section
 - Pins
'''

from math import pi, sqrt
from unum.units import N,mm,m,s,kg,t,Pa
from unum import Unum
kN = Unum.unit('kN' ,1e3*N)
MPa = Unum.unit('MPa',1e6*Pa)
g = 10*m/s/s

def Report(s):
    print(s)

class Material():
    ''' Redefine all variables if to be used 
    with units.'''
    def __init__(   self,
                    fy,
                    E = 200000, 
                    v = 0.3, 
                    rho = 7860):
        self.E = E
        self.v = v
        self.rho = rho
        self.fy = fy

class Section():
    '''
    Defines Area and section modulus from 
    different forms and their dimensions
    Hollow and solid forms: Rectangle, 
    Round, Square. Other forms: Roman I, 
    Roman II. Must be double symetric.
    
    Defined properties:
      A = Area, as in: Fx/A
      Iy, Iz, rz, ry = Second moment of Inertia,
        as in: My*rz/Iy
      Qy, Qz, wy, wz = Something for shear,
        as in: Fz*Qy/Iy/wy
        
    Rounds are a class apart
    Squares are generalized into rectangles
    Roman sections are build up from rectangles
              
    Square / Round      /---Do--/
                      t /-/---/ Di
                      
                        +-------+  
                        | +---+ |   
                        | |   | |  
                        | +---+ |
                        +-------+
                   
    Rectangle           /---Dy--/
                        /-/ ty
                        
                        +-------+  /   /
                        | +---+ |  /tz |
                        | |   | |      Dz
                        | +---+ |      |
                        +-------+      /

    Roman I             /---Dy--/
                           /-/ ty
                           
                        +-------+  /  /
                        +--+ +--+  /  |
                           | |     tz Dz
                        +--+ +--+     |
                        +-------+     /

    Roman II            /---Dy--/
                         /--Do-/
                      ty /-/-/ Di
                      
                        +-------+  /  /
                        ++ +-+ ++  /  |
                         | | | |   tz Dz
                        ++ +-+ ++     |
                        +-------+     /
    '''
    def __init__(   self,
                    form,
                    Do=None,
                    Di=None,
                    t=None,
                    Dy=None,
                    Dz=None,
                    ty=None,
                    tz=None):
        # prep
        if form=="square hollow" or form=="round hollow":
            if not t:
                Dy = Do
                dy = Di
            elif not Do:
                Dy = Di+2*t
                dy = Di
            elif not Di:
                Dy = Do
                dy = Do-2*t
            Dz = Dy
            dz = dy
            ry = Do/2
            rz = ry
            wy = Dy-dy
            wz = Dz-dz
        elif form=="square solid" or form=="round solid":
            Dy = Do
            Dz = Do
            dy = 0
            dz = 0
            ry = Do/2
            rz = Do/2
            wy = Dy-dy
            wz = Dz-dz
        elif form=="rectangle":
            Dy = Dy
            Dz = Dz
            if not ty and not tz:
                ty = 0
                tz = 0
                dy = 0
                dz = 0
            else:
                dy = Dy-ty-ty
                dz = Dz-tz-tz
            ry = Dy/2
            rz = Dz/2
            wy = Dy-dy
            wz = Dz-dz
        elif form=="roman I":
            ry = Dy/2
            rz = Dz/2
            wy = ty
            wz = 2*tz # conservative shortcut
        elif form=="roman II":
            if Di:
                Do = Di + 2*ty
            elif Do:
                Di = Do - 2*ty
            else:
                raise "missing dimension"
            ry = Dy/2
            rz = Dz/2
            wy = 2*ty
            wz = 2*tz
        self.Dz = Dz
        self.Dy = Dy
        self.dz = dz
        self.dy = dy
        self.tz = tz
        self.ty = ty
        self.wz = wz
        self.wy = wy
        self.ry = ry
        self.rz = rz
#        calc(form,Dy,dy,ty,wy,Dz,dz,tz,wz)
#    def calc(form,Dy,dy,ty,wy,Dz,dz,tz,wz):
        if form=="square solid" or form=="square hollow" or form== "rectangle":
            A = Dy*Dz-dy*dz
            Iy = 1/12*(  Dy*Dz**3- dy*dz**3 )
            Iz = 1/12*(  Dz*Dy**3- dz*dy**3 )
            Qz = Dz*Dy/2*Dy/4-dz*dy/2*dy/4
            Qy = Dy*Dz/2*Dz/4-dy*dz/2*dz/4
            # Qz = 1/8 *(  Dz*Dy**2- dz*dy**2 )
            # Qy = 1/8 *(  Dy*Dz**2- dy*dz**2 )
        elif form=="round solid" or form=="round hollow":
            A =  pi  *( (Dy/2)**2-(dy/2)**2 )
            Iy = pi/4*( (Dy/2)**4-(dy/2)**4 )
            Iz = pi/4*( (Dz/2)**4-(dz/2)**4 )
            Qz =  2/3*( (Dy/2)**3-(dy/2)**3 )
            Qy =  2/3*( (Dz/2)**3-(dz/2)**3 )
        elif form=="roman I":
            so = Section("rectangle",Dy=Dy, Dz=Dz)
            sm = Section("rectangle",Dy=Dy, Dz=Dz-tz-tz)
            si = Section("rectangle",Dy=ty, Dz=Dz-tz-tz)
            A  = so.A  - sm.A  + si.A
            Iy = so.Iy - sm.Iy + si.Iy
            Iz = so.Iz - sm.Iz + si.Iz
            Qy = so.Qy - sm.Qy + si.Qy
            Qz = so.Qz - sm.Qz + si.Qz
        elif form=="roman II":
            so = Section("roman I",Dy=Dy,Dz=Dz,ty=Do,tz=tz)
            si = Section("rectangle",Dy=Di,Dz=Dz-tz-tz)
            A  = so.A  - si.A
            Iy = so.Iy - si.Iy
            Iz = so.Iz - si.Iz
            Qy = so.Qy - si.Qy
            Qz = so.Qz - si.Qz
#        file(form,A,Iy,Iz,Qy,Qz,ry,rz,wy,wz)
#    def file(form,A,Iy,Iz,Qy,Qz,ry,rz,wy,wz):
        self.form = form
        self.A = A
        self.Iy = Iy
        self.Iz = Iz
        self.Qy = Qy
        self.Qz = Qz
        assert Iz*wz/Qz <= A, "{}Iz{},wz{},Qz{},A{}".format(form,Iz,wz,Qz,A)
        assert Iy*wy/Qy <= A, "{}Iy{},wy{},Qy{},A{}".format(form,Iy,wy,Qy,A)
    def report(self):
        if self.form=="rectangle":
            Report("Rectangle HxhxBxb: {}x{}x{}x{}".format(self.Dz,self.dz,self.Dy,self.dy))
        elif self.form=="round hollow":
            Report("Tube Dxt: {}x{}".format(self.Dz,self.tz))


class Load():
    '''
    Defines a load vector at a distance vector
    from the origin of a section
    
         x  _z
         ^  /|
        _|_/__  
       / |/  /  
      /  +--/--> y  
     /_____/  

         x  _z
         ^  /|
         | /  
         |/     
         +----> y  
        _|____  _
       / |   /  /|
      /  V  /  z
     /_____/ |/_
    <--y-->
    
    you know..
    '''
    
    def __init__(self,  Fx=0, Fy=0, Fz=0,
                        dx=0, dy=0, dz=0,
                        Mx=0, My=0, Mz=0):
        self.Fx = Fx
        self.Fy = Fy
        self.Fz = Fz
        self.Mx = Mx + Fz*dy - Fy*dz
        self.My = My + Fx*dz - Fz*dx
        self.Mz = Mz + Fy*dx - Fx*dy
        self.dx = dx
    def report(self):
        txt = 'Load:'
        for l,v in zip("Fx Fy Fz Mx My Mz".split(),[self.Fx,self.Fy,self.Fz,self.Mx,self.My,self.Mz]):
            if v != 0:
                txt += " {} {},".format(l,v)
        Report(txt[:-1])

class Member():
    def __init__(self, mat, sec, load, a=0.34,Le = 2.0):
        self.mat = mat
        self.sec = sec
        self.load = load
        self.a = a
        self.Le = Le
        
        # stability
        if load.dx != 0:
            I = min(sec.Iy, sec.Iz)
            Ncr = pi**2*mat.E*I/(Le*load.dx)**2
            lam = sqrt(sec.A*mat.fy/Ncr)
            phi = 0.5*(1+a*(lam-0.2)+lam**2)
            X = min(1,1/(phi+sqrt(phi**2-lam**2))) 
            self.Rb = mat.fy * sec.A * X
            self.stability = max(0,-load.Fx/self.Rb)
        else:
            self.Rb = False
            self.stability = False

        # strength
        sx = load.Fx/sec.A
        sy = load.My*sec.rz/sec.Iy
        sz = load.Mz*sec.ry/sec.Iz
        ty = load.Fy*sec.Qz/sec.Iz/sec.wz
        tz = load.Fz*sec.Qy/sec.Iy/sec.wy
        sx,sy,sz,ty,tz=map(abs,[sx,sy,sz,ty,tz])
        svm1 = sqrt((sx)**2+3*(ty+tz)**2)
        svm2 = sx+sy+sz
        self.strength = max(svm1,svm2) / mat.fy
        
        # flexibility 
        epsilonx = load.Fx/sec.A/mat.E
        self.displacement = epsilonx*load.dx
        
        self.sec.report()
        self.load.report()
        print( "str:",round(self.strength,2), 
               "sta:",round(self.stability,2),
               "dis:",round(self.displacement,2))

def WeldSection(form,a,t,L):
    '''
    2xFW
    2xPP
    1xPP           /\
          /       /  \
         /\      /    \      
        /| a     |\    \    
       / |  \/   | \    \     __
      /  |  /    |  \    /    /
     /   | /     |   \  /    L
    /____|/      |____\/   _/_
    
         /---t---/   
    
    2xPP
    
    /-a-/--t--/-a-/
    
    +---+     +---+ /
    |   |     |   | |
    |   |     |   | |
    |   |     |   | L
   ...
   
    
    
    '''
    if form=="2xPP":
        return Section("rectangle",Dz=L,dz=0,Dy=a+t+a,dy=t)
    if form=="1xPP":
        pass
    if form=="2xFW":
        pass
 
def Weld(mat,sec,load):
    pass

        
    
class Hole():
    def __init__(self, mat, sec, load):
        '''use section hollow rectangle?'''
        pass
def M355(sec,load, a=0.34,Le = 2.0):
    return Member(Material(355), sec, load, a, Le)
def Rect(**args):
    return Section("rectangle",**args)
