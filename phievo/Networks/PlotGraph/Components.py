import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.path import Path
import numpy as np
from math import *


def Bezier(P0,P1,P2):
    return lambda t: (1-t)**2*P0 + 2*(1-t)*t*P1 + t**2*P2


class Interaction:
    """ In the module Graph, an iteraction between node A and node B stands for at least one edge between those two node.
    It is a mean to keep tracks of all the edges that exist between A and B.
    """
    def __init__(self,node1,node2):
        """Creates an iteraction between node1 and node2. Generally launched when the graph creates an edge between node1 and node2 for the first time.

        Args:
            node1 (Node): Node added first to the graph
            node2 (Node): Node added second to the graph

        Returns:
            :class:`Interaction <phievo.Networks.PlotGraph.Components.Interaction>`:Reference to the interaction
        """

        self.node1 = min(node1,node2)
        self.node2 = max(node1,node2)
        self.isAuto = (node1==node2)
        self.edges = {}
        self.numberEdges = 0
        self.receiveEdge = 0

    def add_edge(self,edge):
        """Add an edge to an the existing interaction

        Args:
            edge (:class:`Edge <phievo.Networks.PlotGraph.Components.Edge>`): edge to be added to the list of edge references
        Returns:
            None
        """

        self.edges[self.numberEdges] = edge
        edge.interaction = self
        self.numberEdges += 1

    def get_patches(self,offsets=(0,0)):
        """Run through the interactions edges to create a Matplotlib patch for each of them

        Args:
            offsets (float,float): Size 2 tuple containing the offset to leave between the edges an the node1 and node2.

        Returns:
            [Matplotlib.Patches]: list of Matplotlib patches
        """

        patches = []
        new_iter = sorted(self.edges.items(), key=lambda a: 1-a[1].receiveEdge)

        for i,edge in new_iter:
            if self.isAuto:
                ## Distribute the arrows symmetrically about the node1-node2 axis
                if self.node1==2:
                    patches.append(edge.get_autoPatch(offsets=offsets,num=i))
                else:
                    patches.append(edge.get_autoPatch(offsets=offsets,num=i))
            else:
                ## Distribute the arrows symmetrically about the node1-node2 axis
                angle = ((1-self.numberEdges)*0.5 + i)*0.4
                if edge.nodeFrom.index==self.node2:
                    angle*=-1
                patches.append(edge.get_patch(offsets=offsets,angle=angle))
        return patches

class Edge:
    """ Directed graph edge between two nodes.
    """
    def __init__(self,nodeFrom,nodeTo,label, **kwargs):
        """ Initialise an edge

        Args:
            nodeFrom (:class:`Node <phievo.Networks.PlotGraph.Components.Node>`): Node at which the edge starts
            nodeTo (:class:`Node <phievo.Networks.PlotGraph.Components.Node>`): Node at which the edge ends
            label (string): Edge label
        """
        self.nodeFrom = nodeFrom
        self.nodeTo = nodeTo
        self.label = label
        self.index = None
        self.center = (0,0)
        self.angle = 0
        self.plot_parameters = dict(kwargs)
        self.angles = []
        self.receiveEdge = 0
        self.interaction = None

    def radius(self,theta):
        return 0
    def record_angle(self,angle):
        self.angles.append(angle)

    def get_vector(self,offsets=(0,0),angle=0):
        """ Generate a starting and ending point of the edge's arrow that accomodates the desired space and between the arrow and the nodes given the node shapes.

        Args:
           offsets (float,float): offsets between the arrow and the two nodes
           angle (float):  If angle is 0, the arrow follows a straigh line between two nodes. Otherwise it is a curved line starting and arriving to the node with two opposite angles with respect to the freeAngle value
        Returns:
            (tuple): tuple containing:
                - start (numpy.array): Start of the arrow
                - end (numpy.array): End of the arrow
        """

        direction =   np.array(self.nodeTo.center) - np.array(self.nodeFrom.center)
        L = np.linalg.norm(direction)
        direction = (direction)/L
        if angle != 0:
            Rot = np.array([[cos(angle)-1,-sin(angle)],[sin(angle),cos(angle)-1]])
            dvec = np.dot(Rot,direction)
            directionFrom = direction + dvec
            directionTo = direction - dvec
        else:
            directionTo = directionFrom = direction

        thetaFrom = np.angle(directionFrom[0] + directionFrom[1]*1j,deg=False)
        thetaTo = np.angle(-directionTo[0] - directionTo[1]*1j, deg=False)
        start = np.array(self.nodeFrom.center) + (offsets[0]+self.nodeFrom.radius(thetaFrom))*directionFrom
        end = np.array(self.nodeTo.center) - (offsets[1]+self.nodeTo.radius(thetaTo))*directionTo
        self.nodeFrom.record_angle(thetaFrom)
        self.nodeTo.record_angle(thetaTo)
        if self.nodeFrom.index<self.nodeTo.index:
            self.angle = np.angle(direction[0] + direction[1]*1j,deg=False)
        else:
            self.angle = -np.angle(direction[0] + direction[1]*1j,deg=False)
        self.set_center(self.compute_center(start,end,-angle))
        return start,end

    def get_vector_auto(self,offsets=(0,0),num=0):
        """ Generate a starting and ending point of the edge's arrow that accomodates the desired space and between the arrow and the node given the node shapes. This is an implementation of get_vector for a edge that start and ends at the same edge.

        Args:
            offsets (float,float): offsets between the arrow and the two nodes
            angle (float):  Here the angle cannot be 0. The arrow is a curved line starting and arriving to the node with two opposite angles with respect to the freeAngle value.
        Returns:
            (tuple): tuple containing:
                - start (numpy.array): Start of the arrow
                - end (numpy.array): End of the arrow
        """
        angle = 0.3 + num*0.1
        node = self.nodeFrom
        freeAngle = node.find_freeAngle()
        thetaFrom = freeAngle-angle
        thetaTo = freeAngle+angle

        start = node.center + node.radius(thetaFrom)*np.array([cos(thetaFrom),sin(thetaFrom)])
        end = node.center + node.radius(thetaTo)*np.array([cos(thetaTo),sin(thetaTo)])

        direct = end-start
        ortho = np.array([direct[1],-direct[0]])

        scaling = 1 + num*0.4
        L = scaling*1.15/np.linalg.norm(ortho)
        l = scaling*0.75/np.linalg.norm(direct)
        scalling = 1+1

        P1 = start + L*ortho - l*direct
        P2 = end + L*ortho + l*direct
        return start,P1,P2,end

    def set_center(self,center):
        self.center = center
    def compute_center(self,A,B,angle):
        M = (A+B)/2
        dX = B-A
        dY = np.array([dX[1],-dX[0]])
        C = M + angle*dY
        C = 0.25*A + 0.5*C + 0.25*B
        return C

    def setReceiveEdge(self):
        self.receiveEdge = 1
        if self.interaction:
            self.interaction.receiveEdge = 1

class Arrow(Edge):
    """ The class arrow is inherited from :class:`Networks.PlotGraph.Components.Edge`. It adds extra fonctionalities to generate Matplolib patches.
    """
    def __init__(self,**kwargs):
        """ Initialisation

        Args:
            kwargs (dic): Optional arguments to configure the matplolib FancyArrowPatch's function.
        Returns:
            None
        """
        super(Arrow, self).__init__( **kwargs)

    def get_patch(self,offsets=(0,0),angle=0.2):
        """ Generates a matplotlib patch for the arrow between two nodes. It takes into account the offset to keep between the ends of the arrow and the nodes given the node shapes.

        Args:
            offsets (float,float): offset between node1 and the start of the arrow and offset between node2 and the end of the arrow
        Returns:
            Matplotlib.Patches
        """
        start,end = self.get_vector(offsets,angle)
        return {"pos":self.center,"text":self.label},patches.FancyArrowPatch(start,end,mutation_scale=20,shrinkA=0,shrinkB=0,connectionstyle='arc3, rad=%f'%-angle,**self.plot_parameters)

    def get_autoPatch(self,offsets=(0,0),num=0):
        """ Generates a matplotlib patch for the arrow between two nodes. It takes into account the offset to keep between the ends of the arrow and the node given the node shape. This is an implementation of get_vector for a edge that start and ends at the same node.

        Args:
            offsets (float,float): offset between node and the start of the arrow and offset between node and the end of the arrow
        Returns:
            Matplotlib.Patches
        """
        start,P1,P2,end = self.get_vector_auto(offsets,num)

        verts = [start,P1,P2,end]
        codes = [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,]
        path = Path(verts, codes)

        param = dict(**self.plot_parameters)
        #arrowstyle = param.pop('arrowstyle')
        return {"pos":self.center,"text":self.label} , patches.FancyArrowPatch(path=path,mutation_scale=20,**param)


class Line(Edge):
    def __init__(self,**kwargs):
        super(Arrow, self).__init__( **kwargs)


## Register "-|" as a new arrow style.
class BarB(patches.ArrowStyle._Bracket):
    """
    An arrow with a bar(|) at the B end. The class is added to matplotlib to allow "-\|" style of arrow.
    """

    def __init__(self,widthB=0.4, angleB=None):
        """
        Initialisation
        Args
            widthB (float): width of the bracket
            lengthB (float):length of the bracket
            angleB (float): angle between the bracket and the line
        """

        super(BarB, self).__init__(None, True,widthB=widthB, lengthB=0, angleB=angleB)

patch = patches.ArrowStyle.register("-|",BarB)



### Nodes
class Node:
    """
    Directed graph node or vertex.
    """
    def __init__(self,label,size,*args, **kwargs):
        """
        Initialisation of the node.

        Args:
            label (str): Label of the node. It identifies the node in the graph. It also severs as the node label when the graph is plotted.
            size (float): Relative node area as compare to the default area.
            kwargs (dict): Dictionary  of options to be used by matplotlib when generating the patch.
        Returns:
            :class:`Networks.PlotGraph.Components.Node`: Reference the the newly created node
        """
        self.label = label
        self.size = size
        self.plot_parameters = dict(kwargs)
        self.center = (0,0)
        self.index = None
        self.angles = [] ## To store the angle of the edges in order to optipistion auto-interaction positioning.

    def set_center(self,pos):
        """ Set the coordinates of the node's center.

        Args:
            pos ( list(float) ): Coordinates of the node's center

        Returns:
            None
        """
        self.center = tuple(pos)

    def record_angle(self,angle):
        """
        Every point on boundary of the Node is refered to by an angle. This function records the postition each time a new edge is drawn. The list of angle is used to choose the optimal position where to add looping edges.

        Args:
            angle (float): Value between 0 and 2π where an new edge arrives or leaves the node.
        Returns:
            None
        """
        self.angles.append(angle)

    def find_freeAngle(self):
        """
        Searches for the best position where to add a new edge to the node. It is used only for looping edges. It tries to increase the angle between the new angle and the already plotted edges.

        Args:
            angle (float): Value between 0 and 2π where an new edge arrives or leaves the node.
        Returns:
            float: the function returns the optimal angle
        """

        angles = np.sort(self.angles)
        if len(angles) == 0:
            return pi/2
        elif len(angles) == 1:
            return -angles[0]
        else:
            interAngles = np.append(angles[1::],angles[0]) - angles
            interAngles[-1] += 2*pi
            maxIndex = np.argmax(interAngles)
            freeAngle = angles[maxIndex] + interAngles[maxIndex]/2
            return freeAngle
    def plot_label(self):
        """
        Write the node label on the plot a the node's center.
        """
        plt.text(x=self.center[0],y=self.center[1],s=self.label,verticalalignment='center', horizontalalignment='center',size="large")


class Circle(Node):
    """
    Circle is inherited from :class:`Networks.PlotGraph.Components.Node` and represents a node with a circular shape (⬤).

    """
    def __init__(self,*args, **kwargs):
        """ See Node
        """
        super(Circle, self).__init__(*args, **kwargs)
        self.R = sqrt(self.size/pi)

    def radius(self,theta):
        """ Every point on the node's boundary is refered to by an angle in rad. Given the shape of the node, compute the radius  of the boundary for a angle.

        .. math::
            \\theta \\rightarrow R

        Args:
            theta (float): Angle a which to compute the distance between the center and the boundary.
        Returns:
            float: corresponding to the radius.
        """
        r = self.R
        return r
    def get_patch(self):
        """
        Draw of a matplotlib patch to be added to the graph plot.

        Returns:
            Matplotlib.Patch
        """
        circle = patches.Circle(self.center,self.R,**self.plot_parameters)
        return circle


class HouseUp(Node):
    """
    :class:`Node <phievo.Networks.PlotGraph.Components.Node>` with a pentagon shape (⬟)
    """

    def __init__(self, *args, **kwargs):
        """ See Node.
        """
        super(HouseUp, self).__init__(*args, **kwargs)
        self.R = sqrt(self.size/(0.5*5*sin(2*pi/5))) ## Default area of 1
    def radius(self,theta):
        """ Every point on the node's boundary is refered to by an angle in rad. Given the shape of the node, compute the radius  of the boundary for a angle.

        .. math::
            \\theta \\rightarrow R \\times \\frac{\\cos \\pi/5}{\\cos((5\\theta + 3\\pi/2)\\%(2pi)/5  - \\pi/5)}

        Args:
            theta (float): Angle a which to compute the distance between the center and the boundary.
        Returns:
            float: corresponding to the radius.
        """
        r=self.R*cos(pi/5)/cos((5*theta + 3*pi/2)%(2*pi)/5  - pi/5);
        return r
    def get_patch(self):
        """
        Draw of a matplotlib patch to be added to the graph plot.

        Returns:
            Matplotlib.Patch
        """
        xx = self.center[0] + self.R*np.cos(np.linspace(0,2*pi-2*pi/5,5)+pi/2)
        yy = self.center[1] + self.R*np.sin(np.linspace(0,2*pi-2*pi/5,5)+pi/2)
        points = [[xx[i],yy[i]] for i in range(5)]
        polygon = patches.Polygon(points,closed=True,**self.plot_parameters)
        return(polygon)



class HouseDown(Node):
    """
    :class:`Node <phievo.Networks.PlotGraph.Components.Node>` with a pentagon shape (⯂).
    """
    def __init__(self, *args, **kwargs):
        """ See Node.
        """
        super(HouseDown, self).__init__(*args, **kwargs)
        self.R = sqrt(self.size/(0.5*5*sin(2*pi/5))) ## Default area of 1
    def radius(self,theta):
        """ Every point on the node's boundary is refered to by an angle in rad. Given the shape of the node, compute the radius  of the boundary for a angle.

        .. math::
            \\theta \\rightarrow R \\times \\frac{\\cos \\pi/5}{\\cos((5\\theta - 3\\pi/2)\\%(2pi)/5  - \\pi/5)}

        Args:
            theta (float): Angle a which to compute the distance between the center and the boundary.
        Returns:
            float: corresponding to the radius.
        """
        r=self.R*cos(pi/5)/cos((5*theta - 3*pi/2)%(2*pi)/5  - pi/5);
        return r
    def get_patch(self):
        """
        Draw of a matplotlib patch to be added to the graph plot.

        Returns:
            Matplotlib.Patch
        """
        xx = self.center[0] + self.R*np.cos(np.linspace(0,2*pi-2*pi/5,5)-pi/2)
        yy = self.center[1] + self.R*np.sin(np.linspace(0,2*pi-2*pi/5,5)-pi/2)
        points = [[xx[i],yy[i]] for i in range(5)]
        polygon = patches.Polygon(points,closed=True,**self.plot_parameters)
        return(polygon)


class TriangleUp(Node):
    """
    :class:`Node <phievo.Networks.PlotGraph.Components.Node>` with a triangle shape (▲).
    """
    def __init__(self,*args, **kwargs):
        super(TriangleUp, self).__init__(*args, **kwargs)
        self.R = sqrt(self.size/(0.5*3*sin(2*pi/3))) ## Default area of 1
    def radius(self,theta):
        """ Every point on the node's boundary is refered to by an angle in rad. Given the shape of the node, compute the radius  of the boundary for a angle.

        .. math::
            \\theta \\rightarrow R \\times \\frac{\\cos \\pi/3}{\\cos((3\\theta - 3\\pi/2)\\%(2pi)/3  - \\pi/3)}

        Args:
            theta (float): Angle a which to compute the distance between the center and the boundary.
        Returns:
            float: corresponding to the radius.
        """
        r=self.R*cos(pi/3)/cos((3*theta - 3*pi/2)%(2*pi)/3  - pi/3);
        return r
    def get_patch(self):
        """
        Draw of a matplotlib patch to be added to the graph plot.

        Returns:
            Matplotlib.Patch
        """
        xx = self.center[0] + self.R*np.cos(np.linspace(0,2*pi-2*pi/3,3)+pi/2)
        yy = self.center[1] + self.R*np.sin(np.linspace(0,2*pi-2*pi/3,3)+pi/2)
        points = [[xx[i],yy[i]] for i in range(3)]
        polygon = patches.Polygon(points,closed=True,**self.plot_parameters)
        return(polygon)

class TriangleDown(Node):
    """
    :class:`Node <phievo.Networks.PlotGraph.Components.Node>` with a triangle shape (▼).
    """
    def __init__(self,*args, **kwargs):
        """ See Node.__init__
        """
        super(TriangleDown, self).__init__(*args, **kwargs)
        self.R = sqrt(self.size/(0.5*3*sin(2*pi/3))) ## Default area of 1

    def radius(self,theta):
        """ Every point on the node's boundary is refered to by an angle in rad. Given the shape of the node, compute the radius  of the boundary for a angle.

        .. math::
            \\theta \\rightarrow R \\times \\frac{\\cos \\pi/3}{\\cos((3\\theta + 3\\pi/2)\\%(2pi)/3  - \\pi/3)}

        Args:
            theta (float): Angle a which to compute the distance between the center and the boundary.
        Returns:
            float: corresponding to the radius.
        """
        r=self.R*cos(pi/3)/cos((3*theta + 3*pi/2)%(2*pi)/3  - pi/3);
        return r
    def get_patch(self):
        """
        Draw of a matplotlib patch to be added to the graph plot.

        Returns:
            Matplotlib.Patch
        """
        xx = self.center[0] + self.R*np.cos(np.linspace(0,2*pi-2*pi/3,3)-pi/2)
        yy = self.center[1] + self.R*np.sin(np.linspace(0,2*pi-2*pi/3,3)-pi/2)
        points = [[xx[i],yy[i]] for i in range(3)]
        polygon = patches.Polygon(points,closed=True,**self.plot_parameters)
        return(polygon)

class Square(Node):
    """
    :class:`Node <phievo.Networks.PlotGraph.Components.Node>` with a square shape (◼).
    """
    def __init__(self,*args, **kwargs):
        """ See Node.__init__
        """
        super(Square, self).__init__(*args, **kwargs)
        self.R = sqrt(self.size/(0.5*4*sin(2*pi/4))) ## Default area of 1

    def radius(self,theta):
        """ Every point on the node's boundary is refered to by an angle in rad. Given the shape of the node, compute the radius  of the boundary for a angle.

        .. math::
            \\theta \\rightarrow R \\times \\frac{\\cos \\pi/4}{\\cos((4\\theta + 2\\pi/2)\\%(2pi)/4  - \\pi/4)}

        Args:
            theta (float): Angle a which to compute the distance between the center and the boundary.
        Returns:
            float: corresponding to the radius.
        """
        r=self.R*cos(pi/4)/cos((4*theta + 2*pi/2)%(2*pi)/4  - pi/4);
        return r
    def get_patch(self):
        """
        Draw of a matplotlib patch to be added to the graph plot.

        Returns:
            Matplotlib.Patch
        """
        xx = self.center[0] + self.R*np.cos(np.linspace(0,2*pi-2*pi/4,4)-pi/4)
        yy = self.center[1] + self.R*np.sin(np.linspace(0,2*pi-2*pi/4,4)-pi/4)
        points = [[xx[i],yy[i]] for i in range(4)]
        polygon = patches.Polygon(points,closed=True,**self.plot_parameters)
        return(polygon)

class RoundedRectangle(Node):
    """
    :class:`Node <phievo.Networks.PlotGraph.Components.Node>` with a RoundedRectangle shape (▢).
    """
    def __init__(self,*args, **kwargs):
        """ See Node.__init__
        """
        super(RoundedRectangle, self).__init__(*args, **kwargs)
        a = sqrt(self.size/(0.5*4*sin(2*pi/4))) ## Default area of 1
        self.a = a
        b = a/2
        self.b = b
        rad = a/4
        self.rad = rad
        self.R = sqrt(a**2 + b**2)
        theta1 = atan(b/(a + rad))
        theta2 = atan((b + rad)/a)
        centerTheta = (theta2+theta1)/2
        vecToAngle = self.R*np.array([cos(centerTheta),sin(centerTheta)])
        self.theta_bound = [theta1,theta2,pi- theta2,pi-theta1,pi+theta1,pi+theta2,2*pi-theta2, 2*pi-theta1]
        self.angle_radius = lambda alpha: np.linalg.norm(vecToAngle + rad*np.array([cos(alpha),sin(alpha)]))
    def radius(self,theta):
        """ Every point on the node's boundary is refered to by an angle in rad. Given the shape of the node, compute the radius  of the boundary for a angle.

        Args:
            theta (float): Angle a which to compute the distance between the center and the boundary.
        Returns:
            float: corresponding to the radius.
        """
        theta = theta%(2*pi)
        if self.theta_bound[0] <= theta <= self.theta_bound[1]:
            alpha = pi/2*(theta-self.theta_bound[0])/(self.theta_bound[1] - self.theta_bound[0])
            r = self.angle_radius(alpha)
        elif self.theta_bound[2] <= theta <= self.theta_bound[3]:
            alpha = pi/2 - pi/2*(theta-self.theta_bound[2])/(self.theta_bound[3] - self.theta_bound[2])

            r = self.angle_radius(alpha)
        elif self.theta_bound[4] <= theta <= self.theta_bound[5]:
            alpha = pi/2*(theta-self.theta_bound[4])/(self.theta_bound[5] - self.theta_bound[4])
            r = self.angle_radius(alpha)
        elif self.theta_bound[6] <= theta <= self.theta_bound[7]:
            alpha = pi/2 - pi/2*(theta-self.theta_bound[6])/(self.theta_bound[7] - self.theta_bound[6])
            r = self.angle_radius(alpha)
        elif 0<=theta<self.theta_bound[0] or self.theta_bound[3]<theta<self.theta_bound[4] or self.theta_bound[7]<theta<2*pi:

            r = (self.a+self.rad)/abs(cos(theta))
        else:
            r = (self.b+self.rad)/abs(sin(theta))

        return r
    def get_patch(self):
        """
        Draw of a matplotlib patch to be added to the graph plot.

        Returns:
            Matplotlib.Patch
        """
        rect = patches.FancyBboxPatch(self.center-np.array([self.a,self.b]),2*self.a,2*self.b,boxstyle='round, pad=%f'%self.rad,**self.plot_parameters)
        return(rect)
