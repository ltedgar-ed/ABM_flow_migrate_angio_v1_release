#------------------------------------------------------------------------------
# A module containing classes for various tensors in three-dimensional space
#
# Lowell Taylor Edgar
# Usher Institute of Population Health Sciences and Informatics
# University of Edinburgh
# 2019

from math import sqrt
import numpy

#------------------------------------------------------------------------------
# Vector class in three-dimensions (first-order tensor)
class Vect:
    ###
    # Vector constructor
    def __init__(self, x=0., y=0., z=0.):
        assert (type(x) is int) and (type(y) is int) and (type(z) is int) or (type(x) is float) and (type(y) is float) and (type(z) is float), "Type Error: Vectors can only store int or float"
        
        x : float
        y : float
        z : float
        
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
    
    
    ###
    # Vector operators
    
    # Add a number or another vector
    def __add__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Vect)), "Type Error: Vector addition requires type int, float, or Vect"
        
        if (type(other) is int) or (type(other) is float):
            return Vect(self.x + other, self.y + other, self.z + other)
        
        if (type(other) is Vect):
            return Vect(self.x + other.x, self.y + other.y, self.z + other.z)
    
    # Add to number or another vector
    def __radd__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Vect)), "Type Error: Vector addition requires type int, float, or Vect"
        
        if (type(other) is int) or (type(other) is float):
            return Vect(other + self.x, other + self.y, other + self.z)
        
        if (type(other) is Vect):
            return Vect(other.x + self.x, other.y + self.y, other.z + self.z)
        
    
    # Subtract a number or another vector
    def __sub__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Vect)), "Type Error: Vector subtraction requires type int, float, or Vect"
        
        if (type(other) is int) or (type(other) is float):
            return Vect(self.x - other, self.y - other, self.z - other)
        
        if (type(other) is Vect):
            return Vect(self.x - other.x, self.y - other.y, self.z - other.z)
        
    # Subtract from a number or another vector
    def __rsub__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Vect)), "Type Error: Vector subtraction requires type int, float, or Vect"
        
        if (type(other) is int) or (type(other) is float):
            return Vect(other - self.x, other - self.y, other - self.z)
        
        if (type(other) is Vect):
            return Vect(other.x - self.x, other.y - self.y, other.z - self.z)
    
    
    # Multiply by a number or calculate the dot product with another vector
    def __mul__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Vect)), "Type Error: Vector multiplication requires type int, float, or Vect"
        
        if (type(other) is int) or (type(other) is float):
            return Vect(self.x * other, self.y * other, self.z * other)
        
        # Calculate the dot product if multiplying with another vector
        if (type(other) is Vect):
            return self.x*other.x + self.y*other.y + self.z*other.z
        
    # Multiply by a number or calculate the dot product with another vector
    def __rmul__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Vect)), "Type Error: Vector multiplication requires type int, float, or Vect"
        
        if (type(other) is int) or (type(other) is float):
            return Vect(other * self.x, other * self.y, other * self.z)
        
        # Calculate the dot product if multiplying with another vector
        if (type(other) is Vect):
            return other.x*self.x + other.y*self.y + other.z*self.z
    
    
    # Divide by a number
    def __truediv__(self, other):
        assert ((type(other) is int) or (type(other) is float)), "Type Error: Vector division requires type int or float"
        
        return Vect(self.x / other, self.y / other, self.z / other)
        
    
    # Calculate the negative vector
    def __neg__(self):
        return Vect(-self.x, -self.y, -self.z)
    
    
    # Determine if vector is equal to another vector
    def __eq__(self, other):
        assert (type(other) is Vect), "Type Error: Can only check for equality with another Vect"
        
        if (self.x == other.x) and (self.y == other.y) and (self.z == other.z):
            return True
        else:
            return False 

    
    # Determine if vector is not equal to another vector
    def __ne__(self, other):
        assert (type(other) is Vect), "Type Error: Can only check for equality with another Vect"
        
        if (self.x != other.x) and (self.y != other.y) and (self.z != other.z):
            return True
        else:
            return False
        
        
    # Iterator for the components of the vector
    def __iter__(self):
        vect_iter = iter([self.x, self.y, self.z])
        
        return vect_iter
    
    
    # String conversion function for printing
    def __str__(self):
        return "({0.x}, {0.y}, {0.z})".format(self)
    
    
    # Round the components of the vector
    def __round__(self, n):
        ret = Vect()
        
        ret.x = round(self.x, n)
        ret.y = round(self.y, n)
        ret.z = round(self.z, n)
        
        return ret
    
    
    ###
    # Vector methods
    
    # Create a copy of the vector
    def copy(self):
        return Vect(self.x, self.y, self.z)
    
    
    # Calculate vector length
    def length(self):
        return sqrt(self * self)
    
    
    # Obtain the unit vector
    def unit(self):
        length = self.length()
        
        if (length != 0):
            self.x /= length
            self.y /= length
            self.z /= length
    
    
    # Calculate the outer (tensor) product with another vector
    def outer(self, other):
        assert (type(other) is Vect), "Type Error: Outer product requires two vectors"

        ret = Tensor2O()
        
        ret.xx = self.x * other.x
        ret.xy = self.x * other.y
        ret.xz = self.x * other.z
        ret.yx = self.y * other.x
        ret.yy = self.y * other.y
        ret.yz = self.y * other.z
        ret.zx = self.z * other.x
        ret.zy = self.z * other.y
        ret.zz = self.z * other.z
        
        return ret
        
    
    # Calculate the cross product with another vector
    def cross(self, other):
        assert (type(other) is Vect), "Type Error: Cross product requires another vector"
        
        ret = Vect()
        ret.x = self.y*other.z - self.z*other.y
        ret.y = self.z*other.x - self.x*other.z
        ret.z = self.x*other.y - self.y*other.x
        
        return ret
    
    
    # Apply a change of basis to the vector
    def change_basis(self, exp, eyp, ezp):
        assert ((type(exp) is Vect) and (type(eyp) is Vect) and (type(ezp) is Vect)), "Type Error: Basis change requires 3 new orthogonal unit vectors that make up the new basis"
        
        ex = Vect(1, 0, 0)
        ey = Vect(0, 1, 0)
        ez = Vect(0, 0, 1)
        
        Q = Tensor2O()
        Q.xx = ex * exp
        Q.xy = ex * eyp
        Q.xz = ex * ezp
        
        Q.yx = ey * exp
        Q.yy = ey * eyp
        Q.yz = ey * ezp
        
        Q.zx = ez * exp
        Q.zy = ez * eyp
        Q.zz = ez * ezp
        
        assert (Q.isortho() == True), "Basis Error: Rotation matrix is not orthogonal"
        
        new_vect = Vect()
        new_vect = Q.transpose() * self
        
        self.x = new_vect.x
        self.y = new_vect.y
        self.z = new_vect.z
           

#------------------------------------------------------------------------------
# Tensor class in three-dimensions (second-order tensor)
class Tensor2O:
    ###
    # Tensor constructor
    def __init__(self, components=[0.]*9):
        assert (len(components) == 9), "Constructor Error: Must give 9 components xx, xy, xz, yx, yy, yz, zx, zy, zz"
        
        for comp in components:
            assert ((type(comp) == int) or (type(comp) == float)), "Type Error: Tensor can only store int of float"
        
        xx : float
        xy : float
        xz : float
        yx : float
        yy : float
        yz : float
        zx : float
        zy : float
        zz : float
                
        self.xx = float(components[0])
        self.xy = float(components[1])
        self.xz = float(components[2])
        
        self.yx = float(components[3])
        self.yy = float(components[4])
        self.yz = float(components[5])
        
        self.zx = float(components[6])
        self.zy = float(components[7])
        self.zz = float(components[8])
    
    
    ###
    # Tensor operators
    
    # Right Addition operator
    def __add__(self, other):
        # Add a number to the tensor
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Tensor2O)), "Type Error: Tensor addition requires int, float, or Tensor2O"
        
        if (type(other) is int) or (type(other) is float):
            ret = [0.]*9
            
            ret[0] = self.xx + other
            ret[1] = self.xy + other
            ret[2] = self.xz + other
            ret[3] = self.yx + other
            ret[4] = self.yy + other
            ret[5] = self.yz + other
            ret[6] = self.zx + other
            ret[7] = self.zy + other
            ret[8] = self.zz + other
            
            return Tensor2O(ret)
        
        # Add a tensor to another tensor
        if (type(other) is Tensor2O):
            ret = [0.]*9
            
            ret[0] = self.xx + other.xx
            ret[1] = self.xy + other.xy
            ret[2] = self.xz + other.xz
            ret[3] = self.yx + other.yx
            ret[4] = self.yy + other.yy
            ret[5] = self.yz + other.yz
            ret[6] = self.zx + other.zx
            ret[7] = self.zy + other.zy
            ret[8] = self.zz + other.zz
            
            return Tensor2O(ret)
        
    
    # Left Addition operator
    def __radd__(self, other):
        # Add a number to the tensor
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Tensor2O)), "Type Error: Tensor addition requires int, float, or Tensor2O"
        
        if (type(other) is int) or (type(other) is float):
            ret = [0.]*9
            
            ret[0] = other + self.xx
            ret[1] = other + self.xy
            ret[2] = other + self.xz
            ret[3] = other + self.yx
            ret[4] = other + self.yy
            ret[5] = other + self.yz
            ret[6] = other + self.zx
            ret[7] = other + self.zy
            ret[8] = other + self.zz
            
            return Tensor2O(ret)
        
        # Add a tensor to another tensor
        if (type(other) is Tensor2O):
            ret = [0.]*9
            
            ret[0] = other.xx + self.xx
            ret[1] = other.xy + self.xy
            ret[2] = other.xz + self.xz
            ret[3] = other.yx + self.yx
            ret[4] = other.yy + self.yy
            ret[5] = other.yz + self.yz
            ret[6] = other.zx + self.zx
            ret[7] = other.zy + self.zy
            ret[8] = other.zz + self.zz
            
            return Tensor2O(ret)
    
    
    # Right Subtraction operator
    def __sub__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Tensor2O)), "Type Error: Tensor subtraction requires int, float, or Tensor2O"
        
        return self + -other
    
    
    # Left Subtraction operator
    def __rsub__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Tensor2O)), "Type Error: Tensor subtraction requires int, float, or Tensor2O"
        
        return other + -self
    
    
    # Right Multiplication operator
    def __mul__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Vect) or (type(other) is Tensor2O)), "Type Error: Tensor multiplication requires int, float, Vect, or Tensor2O"
        
        # Multiply by a number
        if (type(other) is int) or (type(other) is float):
            ret = [0.]*9
            
            ret[0] = self.xx * other
            ret[1] = self.xy * other
            ret[2] = self.xz * other
            ret[3] = self.yx * other
            ret[4] = self.yy * other
            ret[5] = self.yz * other
            ret[6] = self.zx * other
            ret[7] = self.zy * other
            ret[8] = self.zz * other
            
            return Tensor2O(ret)
        
        # Tensor-vector multiplication
        if (type(other) is Vect):
            ret = Vect()
            
            ret.x = self.xx*other.x + self.xy*other.y + self.xz*other.z
            ret.y = self.yx*other.x + self.yy*other.y + self.yz*other.z
            ret.z = self.zx*other.x + self.zy*other.y + self.zz*other.z
            
            return ret
        
        # Tensor-tensor multiplication
        if (type(other is Tensor2O)):
            ret = Tensor2O()
            
            ret.xx = self.xx*other.xx + self.xy*other.yx + self.xz*other.zx
            ret.xy = self.xx*other.xy + self.xy*other.yy + self.xz*other.zy
            ret.xz = self.xx*other.xz + self.xy*other.yz + self.xz*other.zz
            
            ret.yx = self.yx*other.xx + self.yy*other.yx + self.yz*other.zx
            ret.yy = self.yx*other.xy + self.yy*other.yy + self.yz*other.zy
            ret.yz = self.yx*other.xz + self.yy*other.yz + self.yz*other.zz
            
            ret.zx = self.zx*other.xx + self.zy*other.yx + self.zz*other.zx
            ret.zy = self.zx*other.xy + self.zy*other.yy + self.zz*other.zy
            ret.zz = self.zx*other.xz + self.zy*other.yz + self.zz*other.zz
            
            return ret
        
    
    # Left Multiplication operator
    def __rmul__(self, other):
        assert ((type(other) is int) or (type(other) is float) or (type(other) is Vect) or (type(other) is Tensor2O)), "Type Error: Tensor multiplication requires int, float, Vect, or Tensor2O"
        
        # Multiply by a number
        if (type(other) is int) or (type(other) is float):
            ret = [0.]*9
            
            ret[0] = other * self.xx
            ret[1] = other * self.xy
            ret[2] = other * self.xz
            ret[3] = other * self.yx
            ret[4] = other * self.yy
            ret[5] = other * self.yz
            ret[6] = other * self.zx
            ret[7] = other * self.zy
            ret[8] = other * self.zz
            
            return Tensor2O(ret)
        
        # Tensor-vector multiplication
        if (type(other) is Vect):
            ret = Vect()
            
            ret.x = other.x*self.xx + other.y*self.yx + other.z*self.zx
            ret.y = other.x*self.xy + other.y*self.yy + other.z*self.zy
            ret.z = other.x*self.xz + other.y*self.yz + other.z*self.zz
            
            return ret
        
        # Tensor-tensor multiplication
        if (type(other is Tensor2O)):
            ret = Tensor2O()
            
            ret.xx = other.xx*self.xx + other.xy*self.yx + other.xz*self.zx
            ret.xy = other.xx*self.xy + other.xy*self.yy + other.xz*self.zy
            ret.xz = other.xx*self.xz + other.xy*self.yz + other.xz*self.zz
            
            ret.yx = other.yx*self.xx + other.yy*self.yx + other.yz*self.zx
            ret.yy = other.yx*self.xy + other.yy*self.yy + other.yz*self.zy
            ret.yz = other.yx*self.xz + other.yy*self.yz + other.yz*self.zz
            
            ret.zx = other.zx*self.xx + other.zy*self.yx + other.zz*self.zx
            ret.zy = other.zx*self.xy + other.zy*self.yy + other.zz*self.zy
            ret.zz = other.zx*self.xz + other.zy*self.yz + other.zz*self.zz
            
            return ret
        
    
    # Division operator
    def __truediv__(self, other):
        assert (type(other) is int) or (type(other) is float), "Type Error: Tensor division requires int or float"

        ret = [0.]*9
        
        ret[0] = self.xx / other
        ret[1] = self.xy / other
        ret[2] = self.xz / other
        ret[3] = self.yx / other
        ret[4] = self.yy / other
        ret[5] = self.yz / other
        ret[6] = self.zx / other
        ret[7] = self.zy / other
        ret[8] = self.zz / other
        
        return Tensor2O(ret)
    
    
    # Negative operator
    def __neg__(self):
        return Tensor2O([-self.xx, -self.xy, -self.xz, -self.yx, -self.yy, -self.yz, -self.zx, -self.zy, -self.zz])
    
    
    # Determine if tensor is equal to another tensor
    def __eq__(self, other):
        assert (type(other) is Tensor2O), "Type Error: Can only check for equality with another Tensor"
        
        if (self.xx == other.xx) and (self.xy == other.xy) and (self.xz == other.xz) and (self.yx == other.yx) and (self.yy == other.yy) and (self.yz == other.yz) and (self.zx == other.zx) and (self.zy == other.zy) and (self.zz == other.zz):
            return True
        else:
            return False
        
    
    # Determine if tensor is equal to another tensor
    def __ne__(self, other):
        assert (type(other) is Tensor2O), "Type Error: Can only check for equality with another Tensor"
        
        if (self.xx != other.xx) and (self.xy != other.xy) and (self.xz != other.xz) and (self.yx != other.yx) and (self.yy != other.yy) and (self.yz != other.yz) and (self.zx != other.zx) and (self.zy != other.zy) and (self.zz != other.zz):
            return True
        else:
            return False
        
        
    # Iterator for the components of the tensor
    def __iter__(self):
        tens_iter = iter([self.xx, self.xy, self.xz, self.yx, self.yy, self.yz, self.zx, self.zy, self.zz])
        
        return tens_iter
    
    
    # String conversion function for printing   
    def __str__(self):
        return "[{0.xx}, {0.xy}, {0.xz}\n {0.yx}, {0.yy}, {0.yz}\n {0.zx}, {0.zy}, {0.zz}]".format(self)
    
    
    # Round the components of the tensor
    def __round__(self, n):
        ret = Tensor2O()
        
        ret.xx = round(self.xx, n)
        ret.xy = round(self.xy, n)
        ret.xz = round(self.xz, n)
        ret.yx = round(self.yx, n)
        ret.yy = round(self.yy, n)
        ret.yz = round(self.yz, n)
        ret.zx = round(self.zx, n)
        ret.zy = round(self.zy, n)
        ret.zz = round(self.zz, n)
        
        return ret
    
    
    ###
    # Tensor methods
    
    # Create a copy of the tensor
    def copy(self):
        return Tensor2O([self.xx, self.xy, self.xz, self.yx, self.yy, self.yz, self.zx, self.zy, self.zz])
    
    
    # Calculate the trace of the tensor
    def trace(self):
        return (self.xx + self.yy + self.zz)
    
    
    # Transpose the tensor
    def transpose(self):
        ret = Tensor2O()
        
        ret.xx = self.xx
        ret.yy = self.yy
        ret.zz = self.zz
        
        ret.xy = self.yx
        ret.yx = self.xy
        ret.xz = self.zx
        ret.zx = self.xz
        ret.yz = self.zy
        ret.zy = self.yz
        
        return ret
    
    
    # Return the symmetric portion of the tensor
    def sym(self):
        return (self + self.transpose())/2.
        
    
    # Return the skew-symmetric portion of the tensor
    def skew(self):
        return (self - self.transpose())/2.
    
    
    # Create an identity tensor
    def eye(self):
        return Tensor2O([1, 0, 0, 0, 1, 0, 0, 0, 1])
  
    
    # Calculate the inner product with another tensor
    def inner(self, other):
        assert (type(other) is Tensor2O), "Type Error: Inner product requires another Tensor"
        
        return (self.xx*other.xx + self.xy*other.xy + self.xz*other.xz + self.yx*other.yx + self.yy*other.yy + self.yz*other.yz + self.zx*other.zx + self.zy*other.zy + self.zz*other.zz)
    
    
    # Calculate the outer product with another tensor
    def outer(self, other):
        assert (type(other) is Tensor2O), "Type Error: Outer product requires another tensor"
        
        return (self.xx*other.xx + self.yx*other.xy + self.zx*other.xz + self.xy*other.yx + self.yy*other.yy + self.zy*other.yz + self.xz*other.zx + self.yz*other.zy + self.zz*other.zz)        
        
    
    # Calculate the magnitude of the tensor
    def mag(self):
        return sqrt(self.inner(self))
    

    # Calculate the volumetric component of the tensor
    def volumetric(self):
        I = Tensor2O().eye()
        
        return (I * (self.trace()/3.))
    
    
    # Calculate the deviatoric component of the tensor
    def deviatoric(self):
        
        return self - self.volumetric()
     
    
    # Determine if the tensor is symmetric
    def issym(self):
        return (self == self.transpose())
    
    
    # Determine if the tensor is orthogonal
    def isortho(self):
        I = Tensor2O().eye()
        
        return (self * self.transpose() == I)
    
    
    # Calculate the determinant of the matrix
    def det(self):
        return (self.xx*self.yy*self.zz + self.xy*self.yz*self.zx + self.xz*self.yx*self.zy - self.xz*self.yy*self.zx - self.xy*self.yx*self.zz - self.xx*self.yz*self.zy)
    
    
    # Calculate the inverse of the tensor
    def inv(self):
        assert (self.det() != 0.), "Matrix is not invertible"
        
        inv = Tensor2O()
        
        inv.xx = MinorMatrix([self.yy, self.yz, self.zy, self.zz]).det()
        inv.xy = MinorMatrix([self.xz, self.xy, self.zz, self.zy]).det()
        inv.xz = MinorMatrix([self.xy, self.xz, self.yy, self.yz]).det()
        
        inv.yx = MinorMatrix([self.yz, self.yx, self.zz, self.zx]).det()
        inv.yy = MinorMatrix([self.xx, self.xz, self.zx, self.zz]).det()
        inv.yz = MinorMatrix([self.xz, self.xx, self.yz, self.yx]).det()
      
        inv.zx = MinorMatrix([self.yx, self.yy, self.zx, self.zy]).det()
        inv.zy = MinorMatrix([self.xy, self.xx, self.zy, self.zx]).det()
        inv.zz = MinorMatrix([self.xx, self.xy, self.yx, self.yy]).det()
        
        inv = inv * (1./self.det()) 
        
        return inv
    
    
    # Calcute the adjugate matrix
    def adjugate(self):
        return self.det() * self.inv()
    
    
    # Calculate the cofactor matrix
    def cofactor(self):
        return self.adjugate().transpose()
    
    
    # Calculate the 3 principal invariants (I invariants) of the tensor
    def prin_invariants(self):
        
        I1 = self.trace()
        I2 = (self.trace()**2 - (self * self).trace())/2.
        I3 = self.det()
        
        return I1, I2, I3
    
    
    # Calculate the 3 main invariants of the tensor
    def main_invariants(self):
        I1, I2, I3 = self.prin_invariants()
        
        J1 = I1
        J2 = I1**2 - 2*I2
        J3 = I1**3 - 3*I1*I2 + 3*I3
        
        return J1, J2, J3
    
    
    # Calculate the eigenvalues and eigenvectors of the tensor, returns a list of eigenvalues and a list of eigenvectors
    def eigen(self):
        eigen_values = [0, 0, 0]
        eigen_vectors = []
        
        numpy_matrix = numpy.array([[self.xx, self.xy, self.xz], [self.yx, self.yy, self.yz], [self.zx, self.zy, self.zz]])
        
        eigen_values, vects = numpy.linalg.eig(numpy_matrix)
        
        v1 = Vect(float(vects[0][0]), float(vects[1][0]), float(vects[2][0]))
        v2 = Vect(float(vects[0][1]), float(vects[1][1]), float(vects[2][1]))
        v3 = Vect(float(vects[0][2]), float(vects[1][2]), float(vects[2][2]))
        
        eigen_vectors.append(v1)
        eigen_vectors.append(v2)
        eigen_vectors.append(v3)
        
        return eigen_values, eigen_vectors
    
    
    # Apply a change of basis to the tensor
    def change_basis(self, exp, eyp, ezp):
        assert ((type(exp) is Vect) and (type(eyp) is Vect) and (type(ezp) is Vect)), "Type Error: Basis change requires 3 new orthogonal unit vectors that make up the new basis"
        
        ex = Vect(1, 0, 0)
        ey = Vect(0, 1, 0)
        ez = Vect(0, 0, 1)
        
        Q = Tensor2O()
        Q.xx = ex * exp
        Q.xy = ex * eyp
        Q.xz = ex * ezp
        
        Q.yx = ey * exp
        Q.yy = ey * eyp
        Q.yz = ey * ezp
        
        Q.zx = ez * exp
        Q.zy = ez * eyp
        Q.zz = ez * ezp
        
        assert (Q.isortho() == True), "Basis Error: Rotation matrix is not orthogonal"
    
        new_tens = Tensor2O()
        new_tens = Q.transpose() * (self * Q)
        
        self.xx = new_tens.xx
        self.xy = new_tens.xy
        self.xz = new_tens.xz
        
        self.yx = new_tens.yx
        self.yy = new_tens.yy
        self.yz = new_tens.yz
        
        self.zx = new_tens.zx
        self.zy = new_tens.zy
        self.zz = new_tens.zz
    

#------------------------------------------------------------------------------
# Minor Matrix class for calculating the inverse of a tensor (basically a matrix in 2D)
class MinorMatrix:
    ###
    # MinorMatrix constructor
    def __init__(self, components=[0]*4):
        assert (len(components) == 4), "Constructor Error: Must give 4 components a, b, c, d"
        
        for comp in components:
            assert ((type(comp) == int) or (type(comp) == float)), "Type Error: Tensor can only store int of float"
        
        a : float
        b : float
        c : float
        d : float
        
        self.a = float(components[0])
        self.b = float(components[1])
        self.c = float(components[2])
        self.d = float(components[3])
     
    
    ###    
    # Minor Matrix methods     
    # Calculate the determinant of the minor matrix
    def det(self):
        return self.a*self.d - self.b*self.c
    

#------------------------------------------------------------------------------
# main program
if __name__ == "__main__":
    a = Vect(1, 1, 1)
    b = Vect(1, 2, 3)
    ex = Vect(1, 0, 0)
    
    print(ex.length())
    
    print(a.length())
    print(b)
    b.unit()
    print(b)
    
    
    c = a + 1
    d = b - 2
    e = b + a
    f = b - a
    
    print()
    print(c)
    print(d)
    print(e)
    print(f)
    
    print()
    print(a * 2)
    print(a * b)
    print(b * ex)
    
    print()
    A = Tensor2O()
    print(A)
    
    B = Tensor2O([1, 2, 3, 4, 5, 6, 7, 8, 9])
    print(B.transpose())
    
    I = Tensor2O().eye()
    print(I)
    
    print(B.sym())
    print(B.skew())
    print(B.sym() + B.skew())
    
    C = Tensor2O([1, 2, 3, 4, 5, 6, 7, 8, 0])
    print(C.inv())
    
    print()
    D = Tensor2O([1, 2, 3, 3, 2, 1, 1, 0, -1])
    w, v = D.eigen()
    print(w)
    
    for vect in v:
        print(vect)
    
    print()
    
    for component in b:
        print(component)
    
    print()
    for component in B:
        print(component)
    
    print()
    ey = Vect(0, 1, 0)
    ez = Vect(0, 0, 1)
    
    ex_old = Vect(1, 0, 0)
    
    ex_old.change_basis(ez, ey, ex*(-1.0))
    print(ex_old)
    
    print()
    B.change_basis(ez, ey, ex*(-1.0))
    print(B)