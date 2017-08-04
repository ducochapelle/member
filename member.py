''' 
Todo:
 - Welds
 - Section: Plastic
 - Pins
 - Section: Torsion
'''

from math import pi # pint doesn't like sqrt

def Report(*s):
    ''' gets a [str,val,val,...,val]. If units are used, some
    formatting is done, otherwise not.'''
    txt = s[0]
    if type(s[1]) not in [type(1),type(1.)]:
        vs = [x.to_base_units() for x in s[1:]]
        print(txt.replace("{}","{:~.4g}").format(*vs))
    else:
        vs = s[1:]
        print(txt.replace("{}","{:.4g}").format(*vs))
              
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
            tz = Dz/2-dz/2
            ty = Dy/2-dy/2
        elif form=="square solid" or form=="round solid":
            Dy = Do
            Dz = Do
            dy = Do-Do
            dz = Do-Do
            ry = Do/2
            rz = Do/2
            wy = Dy-dy
            wz = Dz-dz
        elif form=="rectangle":
            Dy = Dy
            Dz = Dz
            if not ty and not tz:
                ty = Dz-Dz
                tz = Dz-Dz
                dy = Dz-Dz
                dz = Dz-Dz
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
            dy = Dz-Dz
            dz = Dz-Dz
        elif form=="roman II":
            if Di:
                Do = Di + 2*ty
                dy = Do/2 + Di/2
            elif Do:
                Di = Do - 2*ty
                dy = Do/2 + Di/2
            else:
                raise "missing dimension"
            ry = Dy/2
            rz = Dz/2
            wy = 2*ty
            wz = 2*tz
            dz = Dz-Dz
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
            # Jx = 
            # Qz = 1/8 *(  Dz*Dy**2- dz*dy**2 )
            # Qy = 1/8 *(  Dy*Dz**2- dy*dz**2 )
        elif form=="round solid" or form=="round hollow":
            A =  pi  *( (Dy/2)**2-(dy/2)**2 )
            Iy = pi/4*( (Dy/2)**4-(dy/2)**4 )
            Iz = pi/4*( (Dz/2)**4-(dz/2)**4 )
            Qz =  2/3*( (Dy/2)**3-(dy/2)**3 )
            Qy =  2/3*( (Dz/2)**3-(dz/2)**3 )
            # Jx = pi/2*( (Dy/2)**4-(dy/2)**4 )
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
            Report("Rectangle HxhxBxb: {} x {} x {} x {}",self.Dz,self.dz,self.Dy,self.dy)
        elif self.form=="round hollow":
            Report("Tube Dxt: {} x {}",self.Dz,self.tz)
        elif self.form=="round solid":
            Report("Rod D: {}",self.Dz)
        elif self.form=="square solid":
            Report("Bar D: {}",sef.Dz)
        elif self.form=="roman II":
            Report("II Beam: {} x {} x {} x {} x {}",self.Dy,self.ty,self.Dz,self.tz,self.dy)
        else:
            Report("No report")

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
    
    def __init__(self,  Fx=None, Fy=None, Fz=None,
                        dx=None, dy=None, dz=None,
                        Mx=None, My=None, Mz=None):
        force = Fx if Fx else Fy if Fy else Fz if Fz else 0
        force = force - force
        moment = Mx if Mx else My if My else Mz if Mz else 0
        moment = moment - moment
        distance = dx if dx else dy if dy else dz if dz else 0
        distance = distance - distance
        self.Fx = Fx if Fx else force
        self.Fy = Fy if Fy else force
        self.Fz = Fz if Fz else force
        self.Mx = Mx if Mx else moment
        self.My = My if My else moment
        self.Mz = Mz if Mz else moment
        self.dx = dx if dx else distance
        self.dy = dy if dy else distance
        self.dz = dz if dz else distance
    def report(self):
        txt = 'Load:'
        vs = []
        for l,v in zip("Fx Fy Fz Mx My Mz".split(),[self.Fx,self.Fy,self.Fz,self.Mx,self.My,self.Mz]):
            if v != 0:
                txt += " "+l+" {},"
                vs += [v]
        Report([txt[:-1]]+vs)

class Member():
    def __init__(self, mat, sec, load, a=0.34,Le = 2.0):
        self.mat = mat
        self.sec = sec
        self.load = load
        self.a = a
        self.Le = Le
        
        # stability
        if load.dx != 0 and load.Fx != 0:
            I = min(sec.Iy, sec.Iz)
            Ncr = pi**2*mat.E*I/(Le*load.dx)**2
            lam = (sec.A*mat.fy/Ncr)**.5
            phi = 0.5*(1+a*(lam-0.2)+lam**2)
            X = min(1,1/(phi+(phi**2-lam**2)**.5)) 
            self.Rb = mat.fy * sec.A * X
            self.stability = max(0,-load.Fx/self.Rb)
        else:
            self.Rb = False
            self.stability = 0

        # strength
        if load.dx * load.dy * load.dz != 0:
            Mx = load.Mx + load.Fz*load.dy - load.Fy*load.dz
            My = load.My + load.Fx*load.dz - load.Fz*load.dx
            Mz = load.Mz + load.Fy*load.dx - load.Fx*load.dy
        else:
            Mx,My,Mz = 0,0,0

        sx = load.Fx/sec.A
        sy = (My*sec.rz/sec.Iy) if My != 0 else 0
        sz = (Mz*sec.ry/sec.Iz) if Mz != 0 else 0
        ty = load.Fy*sec.Qz/sec.Iz/sec.wz
        tz = load.Fz*sec.Qy/sec.Iy/sec.wy
        sx,sy,sz,ty,tz=map(abs,[sx,sy,sz,ty,tz])
        svm1 = ((sx)**2+3*(ty+tz)**2)**.5
        svm2 = sx+sy+sz
        self.strength = max(svm1,svm2) / mat.fy
        
        # flexibility 
        epsilonx = load.Fx/sec.A/mat.E
        self.displacement = epsilonx

        self.stability *= self.strength-self.strength
        self.sec.report()
        # self.load.report()
        Report( "str: {} sta: {} dis: {}",self.strength,self.stability,self.displacement)

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
        return Section("rectangle",Dz=L,Dy=a)
    if form=="2xFW":
        pass        
    
class Hole():
    def __init__(self, mat, sec, load):
        '''use section hollow rectangle?'''
        pass
def M355(sec,load, a=0.34,Le = 2.0):
    return Member(Material(355), sec, load, a, Le)
def Rect(**args):
    return Section("rectangle",**args)

