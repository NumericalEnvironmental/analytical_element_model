##############################################################################################
#
# AEM - semi-analytical solution for steady-state groundwater flow
# in the presence of injectors (point or line sources) and semi-permeable
# fault segments
#
# injectors and faults inherit parameters and methods from generalized linear segment class
#
##############################################################################################

from numpy import *
from scipy.integrate import quad
from scipy.special import *
from scipy.spatial import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from Tkinter import *

##########################################
#
# supporting functions
#
##########################################


dl = 0.01           # offset distance from fault; location at which zero-velocity point will be defined
dx = 0.001          # infinitesimal distance, used to compute gradients


def ReadElements():                                                         ##### read elements definition file and assign objects
    well = []
    injector = []
    fault_int = []
    i = 0
    input_file = open('elements.txt','r')
    for line in input_file:
        if i:                                                               # don't parse header
            line_input = line.split()
            x0 = float(line_input[1])                                       # starting coordinates are in the same column for all element types)
            y0 = float(line_input[2])
            if line_input[0] == 'well':
                Q = float(line_input[3])
                well.append(Well(x0,y0,Q))
            elif line_input[0] == 'injector':
                xf = float(line_input[3])
                yf = float(line_input[4])                
                Qx = float(line_input[5])
                injector.append(HorizInjector(x0,y0,xf,yf,Qx))
            else:               # fault element
                xf = float(line_input[3])
                yf = float(line_input[4])                
                R = float(line_input[5])
                N = int(line_input[6])
                fault_int.append(FaultIntegrated(x0,y0,xf,yf,R,N))
        i += 1
    input_file.close()
    return well,injector,fault_int

def WeightMatrix(well,injector,fault,hydro):                                    ##### construct matrix to solve for strength terms for fault element(s) in model
    num_faults = len(fault)
    A = zeros((num_faults,num_faults),float)
    B = zeros(num_faults,float)
    for row in xrange(num_faults):
        for col in xrange(num_faults): A[row,col] = fault[row].DirGrad(fault[col],hydro)     # impact of fault[col] on fault[row]
        for i in xrange(len(well)): B[row] -= fault[row].DirGrad(well[i],hydro)
        for i in xrange(len(injector)): B[row] -= fault[row].DirGrad(injector[i],hydro)
    return A,B


##########################################
#
# classes
#
##########################################

class Track:
    def __init__(self):         # read in basic attributes for particle tracking routine
        line_input = []     
        input_file = open('track.txt','r')
        for line in input_file: line_input.append(line.split())
        input_file.close()
        self.track = int(line_input[0][1])                      # 1 = do particle tracking, 0 = skip
        self.dir = int(line_input[1][1])                        # 1 = forward tracking, -1 = reverse tracking
        self.t_end = float(line_input[2][1])
        self.t_steps = int(line_input[3][1])
        self.R = float(line_input[4][1])                        # retardation coefficient
        self.t = linspace(0.,self.t_end,self.t_steps)           # output times for particle positions
        self.Interface()
        if self.track:
            self.particle = self.ReadParticleStart()            # particle = n x 2 array of starting particle positions
            output_file = open('particle_tracks.txt','w')       # open output file and write header
            output_file.writelines(['x','\t','y','\t','t','\t','particle','\n'])
            output_file.close()
    def ReadParticleStart(self):
        # read in initial particle positions
        line_input = []     
        input_file = open('particles.txt','r') 
        for line in input_file: line_input.append(line.split())
        input_file.close()
        particle = zeros((len(line_input)-1,2),float)
        for i in xrange(1,len(line_input)):
            for j in xrange(2):
                particle[i-1,j] = float(line_input[i][j])
        return particle
    def Vel(self,part,t,base,well,injector,fault,hydro):
        # components of velocity field at point p, i.e. dx/dt and dy/dt
        x = part[0]
        y = part[1]
        dhdx = (base.Head(x+dx,y,well,injector,fault,hydro) - base.Head(x-dx,y,well,injector,fault,hydro))/(2.*dx)
        dhdy = (base.Head(x,y+dx,well,injector,fault,hydro) - base.Head(x,y-dx,well,injector,fault,hydro))/(2.*dx)
        vx = -hydro.K*dhdx/(hydro.phi*self.R) * sign(self.dir)     
        vy = -hydro.K*dhdy/(hydro.phi*self.R) * sign(self.dir)
        return [vx,vy]
    def MoveODE(self,base,well,injector,fault,hydro):
        # solve transport ODE for particle set and write to output file
        trail = []
        for i,part in enumerate(self.particle):
            position = odeint(self.Vel,part,self.t,args=(base,well,injector,fault,hydro))
            trail.append(transpose(position))                   # record particle trail (to be returned by method for later plotting)
            output_file = open('particle_tracks.txt','a')       # open output file for append
            for j in xrange(len(self.t)): output_file.writelines([str(position[j,0]),'\t',str(position[j,1]),'\t',str(self.t[j]),'\t',str(i),'\n'])
            output_file.close()
        return trail
    def Interface(self):                                                    # create a data input window for modifying hydro parameter values
        root = Tk()
        container = Frame(root)
        container.grid()
        Label(container,text='Particle Tracking',font=('Courier',14)).grid()        
        labels = ['Time span','No. of steps','Retardation coefficient']
        params = [self.t_end,self.t_steps,self.R]
        entry = []
        # check box
        self.w = IntVar()
        self.c_box = Checkbutton(container,text="Track particles",variable=self.w)
        self.c_box.grid(sticky=W)
        self.w.set(self.track)        
        # radio button
        modes = [('Forward',1),('Reverse',-1)]
        self.v = IntVar()
        for item,switch in modes:
            self.r_button = Radiobutton(container,text=item,variable=self.v,value=switch)
            self.r_button.grid(sticky=W)
        self.v.set(self.dir)
        # text boxes
        for i in xrange(len(labels)):
            Label(container,text=labels[i]).grid(row=i+4,sticky=W)
            entry.append(Entry(container))
            entry[i].grid(row=i+4,column=1)
            entry[i].insert(0,str(params[i]))
        btn = Button(container,text='UPDATE',command=lambda:self.Button_click(entry)).grid(column=1)
        root.mainloop()
    def Button_click(self,entry):
        self.track = self.w.get()                      # 1 = do particle tracking, 0 = skip
        self.dir = self.v.get()                        # 1 = forward tracking, -1 = reverse tracking
        self.t_end = float(entry[0].get())
        self.t_steps = int(entry[1].get())
        self.R = float(entry[2].get())        
        self.WriteParams()        
    def WriteParams(self):                                               
        # write present value set to text file
        output_file = open('track.txt','w')
        output_file.writelines(['track','\t',str(self.track),'\n'])
        output_file.writelines(['direction','\t',str(self.dir),'\n'])
        output_file.writelines(['time_span','\t',str(self.t_end),'\n'])
        output_file.writelines(['steps','\t',str(self.t_steps),'\n'])
        output_file.writelines(['R','\t',str(self.R),'\n'])        
        output_file.close()

class Base:
    def __init__(self):                     # read base input file (grid + ambient GW surface)
        line_input = []     
        input_file = open('base.txt','r')
        for line in input_file: line_input.append(line.split())
        input_file.close()
        self.start = [float(line_input[1][1]),float(line_input[1][2])]
        self.end = [float(line_input[2][1]),float(line_input[2][2])]		
        self.N = [int(line_input[3][1]),int(line_input[3][2])]				
        self.x_grid = linspace(self.start[0],self.end[0],self.N[0])  # x-gridding
        self.y_grid = linspace(self.start[1],self.end[1],self.N[1])  # y-gridding
        self.origin = [float(line_input[4][1]),float(line_input[4][2])]
        self.grad_0 = [float(line_input[5][1]),float(line_input[5][2])]
        self.Interface()
    def Interface(self):                    # bring up interface to enable modification
        root = Tk()
        container = Frame(root)
        container.grid()
        Label(container,text='Grid Attributes', font=('Courier',14)).grid()
        column_labels = ['X','Y']
        row_labels = ['Start','End','N','Origin','Gradient']
        for i in xrange(2): Label(container,text=column_labels[i]).grid(column=i+1)
        for i in xrange(5): Label(container,text=row_labels[i]).grid(row=i+1,sticky=W)
        entry_start = []
        entry_end = []
        entry_N = []
        entry_origin = []
        entry_grad_0 = []        
        for i in xrange(2):
            entry_start.append(Entry(container))
            entry_start[i].grid(row=1,column=i+1)
            entry_start[i].insert(0,str(self.start[i]))
        for i in xrange(2):
            entry_end.append(Entry(container))
            entry_end[i].grid(row=2,column=i+1)
            entry_end[i].insert(0,str(self.end[i]))
        for i in xrange(2):
            entry_N.append(Entry(container))
            entry_N[i].grid(row=3,column=i+1)
            entry_N[i].insert(0,str(self.N[i]))
        for i in xrange(2):
            entry_origin.append(Entry(container))
            entry_origin[i].grid(row=4,column=i+1)
            entry_origin[i].insert(0,str(self.origin[i]))
        for i in xrange(2):
            entry_grad_0.append(Entry(container))
            entry_grad_0[i].grid(row=5,column=i+1)
            entry_grad_0[i].insert(0,str(self.grad_0[i]))
        btn = Button(container,text='UPDATE',command=lambda:self.Button_click(entry_start,entry_end,entry_N,entry_origin,entry_grad_0)).grid(column=2)
        root.mainloop()
    def Button_click(self,entry_start,entry_end,entry_N,entry_aniso,entry_slope):
        self.start = array([float(entry_start[0].get()),float(entry_start[1].get())])
        self.end = array([float(entry_end[0].get()),float(entry_end[1].get())])
        self.N = array([int(entry_N[0].get()),int(entry_N[1].get())])
        self.origin = array([int(entry_origin[0].get()),int(entry_origin[1].get())])
        self.grad_0 = array([int(entry_grad_0[0].get()),int(entry_grad_0[1].get())])
        self.WriteParams()
    def WriteParams(self):
        # write present value set to text file
        output_file = open('base.txt','w')
        output_file.writelines(['\t','X','\t','Y','\n'])
        output_file.writelines(['start','\t',str(self.start[0]),'\t',str(self.start[1]),'\n'])
        output_file.writelines(['end','\t',str(self.end[0]),'\t',str(self.end[1]),'\n'])
        output_file.writelines(['N','\t',str(self.N[0]),'\t',str(self.N[1]),'\n'])
        output_file.writelines(['origin','\t',str(self.origin[0]),'\t',str(self.origin[1]),'\n'])
        output_file.writelines(['gradient','\t',str(self.grad_0[0]),'\t',str(self.grad_0[1]),'\n'])
        output_file.close()
    def Head(self,x,y,well,injector,fault,hydro):
        # superimpose head-change components
        h = hydro.h0 + self.grad_0[0]*(x-self.origin[0]) + self.grad_0[1]*(y-self.origin[1])
        for i in xrange(len(well)): h += well[i].DeltaH(x,y,hydro)
        for i in xrange(len(injector)): h += injector[i].DeltaH(x,y,hydro)
        for i in xrange(len(fault)): h += fault[i].DeltaH(x,y,hydro)
        return h
    def ProcessOutput(self,well,injector,fault,hydro,trail):
        h = zeros(len(self.x_grid)*len(self.y_grid),float)
        # write output file
        output_file = open('heads.txt','w')
        line_out = ['x','\t','y','\t','h','\n']
        output_file.writelines(line_out)
        k = 0
        for y in self.y_grid:
            for x in self.x_grid:
                h[k] = self.Head(x,y,well,injector,fault,hydro)
                line_out = [str(x),'\t',str(y),'\t',str(h[k]),'\n']
                output_file.writelines(line_out)
                k += 1
        output_file.close()
        # create contour visualization
        p = h.reshape(len(self.x_grid),len(self.y_grid))
        x_label = 'x'
        y_label = 'y'
        plt.xlim([self.start[0], self.end[0]])
        plt.ylim([self.start[1], self.end[1]])
        plt.pcolor(self.x_grid,self.y_grid,p,cmap=cm.gray)
        plt.contour(self.x_grid,self.y_grid,p,20)      
        plt.colorbar()
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        # draw particle tracks
        for point_set in trail: plt.plot(point_set[0],point_set[1],'ro')
        plt.show()

class ElementInterface:
    def __init__(self,base,well,injector,fault_int):
        # basic attributes
        root = Tk()                                                                         # container for element attributes
        self.x_span = base.end[0] - base.start[0]
        self.y_span = base.end[1] - base.start[1]
        self.map_scale = 500                                                                # defined ...
        self.height_f = self.y_span/self.x_span
        self.element_group = well + injector + fault_int                                    # combined list of elements
        self.num_wells = len(well)
        self.num_injectors = len(injector)
        self.num_faults = len(fault_int)
        # create window to manage element properties
        container = Frame(root)
        container.grid()
        Label(container,text='Elements',font=('Courier',14)).grid(column=1)
        # create map of domain to display along with input window
        child_window = Toplevel(root)                                       # container for element map
        domain = Canvas(child_window,width=self.map_scale,height=self.height_f*self.map_scale)
        domain.configure(background='deep sky blue')
        domain.grid()
        # set up element type selection radio buttons
        modes = [('Well',0),('Injector',1),('Fault',2)]
        self.v = IntVar()
        for item,switch in modes:
            self.r_button = Radiobutton(container,text=item,variable=self.v,value=switch,command=self.BlankForm)
            self.r_button.grid(sticky=W)
        # set up element attributes form
        self.label = []
        self.entry = []
        self.text_v = []
        self.attribs = ['x-start','y-start','x-finish','y-finish','value','N']        
        for i in xrange(len(self.attribs)): self.text_v.append(StringVar())
        for i in xrange(len(self.attribs)):
            self.text_v[i].set(self.attribs[i])
            self.label.append(Label(container,textvariable=self.text_v[i]).grid(row=i+4,sticky=W))
            self.entry.append(Entry(container))
            self.entry[i].grid(row=i+4,column=1)
        # populate with first element
        self.pointer = 0        
        self.UpdateElementForm()
        # define action buttons
        back_button = Button(container,text='<<',command=lambda:self.Button_click(0,domain,base)).grid(row=10,column=2,sticky=E)
        forward_button = Button(container,text='>>',command=lambda:self.Button_click(1,domain,base)).grid(row=10,column=3)        
        del_button = Button(container,text='DEL',command=lambda:self.Button_click(2,domain,base)).grid(row=10,column=4)  
        update_button = Button(container,text='UPDATE',command=lambda:self.Button_click(3,domain,base)).grid(row=10,column=5)
        file_button = Button(container,text='SAVE',command=lambda:self.Button_click(4,domain,base)).grid(row=10,column=6) 
        self.PopulateMap(domain,base)
        root.mainloop()
    def Button_click(self,i_button,domain,base):            # action button responses
        if i_button == 0:
            # step backward in element_group list
            if self.pointer > 0: self.pointer -= 1
            self.UpdateElementForm()
        elif i_button == 1:
            # step forward in element group list
            if self.pointer < len(self.element_group)-1: self.pointer += 1
            self.UpdateElementForm()
        elif i_button == 2:
            # remove element from element_group list
            if len(self.element_group) > 1:
                del self.element_group[self.pointer]
                self.pointer = len(self.element_group)-1
                self.UpdateElementForm()
                self.PopulateMap(domain,base)
        elif i_button == 3:
            # update the element_group list with new element (taken from current window), note that elements FILE is NOT automatically updated when this button is pressed
            if self.v.get() == 0: self.element_group.append(Well(float(self.entry[0].get()),float(self.entry[1].get()),float(self.entry[4].get())))                                                                  # well
            elif self.v.get() == 1: self.element_group.append(Well(float(self.entry[0].get()),float(self.entry[1].get()),float(self.entry[2].get()),float(self.entry[3].get()),float(self.entry[4].get())))         # horizontal injector
            else: self.element_group.append(Well(float(self.entry[0].get()),float(self.entry[1].get()),float(self.entry[2].get()),float(self.entry[3].get()),float(self.entry[4].get()),int(self.entry[5].get())))  # fault
            self.PopulateMap(domain,base)
        else:
            # write to the elements file
            output_file = open('elements.txt','w')
            output_file.writelines(['type','\t','x-start','\t','y-start','\t','x-finish','\t','y-finish','\t','value','\t','N','\n'])       # header
            for item in elements_group:
                line_output = []
                line_output.append('\t')
                line_output.append(str(item.x0))
                line_output.append('\t')
                line_output.append(str(item.y0))                
                if isinstance(item,Well):
                    line_output.insert(0,'well')
                    line_output.append('\t')
                    line_output.append(str(item.Q))
                elif isinstance(item,HorizInjector):
                    line_output.insert(0,'injector')
                    line_output.append('\t')
                    line_output.append(str(item.xf))
                    line_output.append('\t')
                    line_output.append(str(item.yf))   
                    line_output.append('\t')
                    line_output.append(str(item.Qx))
                else:                           # integrated fault
                    line_output.insert(0,'fault')
                    line_output.append('\t')
                    line_output.append(str(item.xf))
                    line_output.append('\t')
                    line_output.append(str(item.yf))  
                    line_output.append('\t')
                    line_output.append(str(item.R))
                    line_output.append('\t')
                    line_output.append(str(item.N))
                line_output.append('\n')
                output_file.writelines(line_output)
            output_file.close()
    def BlankForm(self):
        # create a blank element form, with configuration based on radiobutton selection for well, injector, or fault
        for i in xrange(len(self.attribs)): self.entry[i].delete(0,END)
        self.entry[0].insert(0,'0')
        self.entry[1].insert(0,'0')
        self.entry[4].insert(0,'0')                    
        if self.v.get() == 0:
            # well format
            self.entry[2].insert(0,'n/a')
            self.entry[3].insert(0,'n/a')
            self.entry[5].insert(0,'n/a')
            self.text_v[4].set('Q')  
        elif self.v.get() == 1:
            # horizontal injector format
            self.entry[2].insert(0,'0')
            self.entry[3].insert(0,'0')
            self.entry[5].insert(0,'n/a') 
            self.text_v[4].set('Q/length')
        else:
            # integrated fault format
            self.entry[2].insert(0,'0')
            self.entry[3].insert(0,'0')
            self.entry[5].insert(0,'1')             
            self.text_v[4].set('Rel. resistance')
    def UpdateElementForm(self):
        # populate element form with attributes from existing element_group set
        for i in xrange(len(self.attribs)): self.entry[i].delete(0,END)
        self.entry[0].insert(0,self.element_group[self.pointer].x0)
        self.entry[1].insert(0,self.element_group[self.pointer].y0)
        if isinstance(self.element_group[self.pointer],Well):
            self.v.set(0)
            self.entry[2].insert(0,'n/a')
            self.entry[3].insert(0,'n/a')
            self.entry[4].insert(0,self.element_group[self.pointer].Q)            
            self.entry[5].insert(0,'n/a')
            self.text_v[4].set('Q')            
        elif isinstance(self.element_group[self.pointer],HorizInjector):
            self.v.set(1)
            self.entry[2].insert(0,self.element_group[self.pointer].xf)
            self.entry[3].insert(0,self.element_group[self.pointer].yf)
            self.entry[4].insert(0,self.element_group[self.pointer].Qx)             
            self.entry[5].insert(0,'n/a') 
            self.text_v[4].set('Q/length')
        else:                           # integrated fault
            self.v.set(2)
            self.entry[2].insert(0,self.element_group[self.pointer].xf)
            self.entry[3].insert(0,self.element_group[self.pointer].yf)
            self.entry[4].insert(0,self.element_group[self.pointer].R)  
            self.entry[5].insert(0,self.element_group[self.pointer].N) 
            self.text_v[4].set('Rel. resistance')
    def CoordTransform(self,x,y,base):
        # convert absolute coordinates into canvas map coordinates
        x_map = self.map_scale * (x-base.start[0])/self.x_span
        y_map = self.height_f*self.map_scale * (1. - (y-base.start[1])/self.y_span)
        return x_map,y_map
    def PopulateMap(self,domain,base):
        # clear map of any old objects and (re)draw (updated) objects
        domain.delete(ALL)       
        for item in self.element_group:
            x_start,y_start = self.CoordTransform(item.x0,item.y0,base)
            if isinstance(item,Well):
                r = self.map_scale/100.
                domain.create_oval(x_start-r,y_start-r,x_start+r,y_start+r)
            elif isinstance(item,HorizInjector):
                x_end,y_end = self.CoordTransform(item.xf,item.yf,base)                
                domain.create_line(x_start,y_start,x_end,y_end)
            else:           # fault
                x_end,y_end = self.CoordTransform(item.xf,item.yf,base)
                domain.create_line(x_start,y_start,x_end,y_end,dash=(4, 4))

class Hydro:                                                                ##### system hydrogeological parameters and methods
    def __init__(self):
        param = []
        input_file = open('hydro.txt','r')
        for line in input_file:
                line_input = line.split()
                param.append(line_input[1])        
        self.K = float(param[0])                                            # aquifer hydraulic conductivity
        self.b = float(param[1])                                            # aquifer thickness
        self.K_c = float(param[2])                                          # cap rock vertical hydraulic conductivity
        self.b_c = float(param[3])                                          # cap rock thickness
        self.phi = float(param[4])                                          # aquifer porosity
        self.h0 = float(param[5])                                           # reference ambient GW head (at base.origin)        
        input_file.close()
        self.B = sqrt(self.K*self.b*self.b_c/self.K_c)                      # leakance term
        self.Interface()
    def LeakyPointSource(self,xw,yw,x,y,Q):                                 ### return the steady-state drawdown at a point associated with a well in a leaky aquifer (Jacob, 1946 solution)
        r = distance.pdist(array([[xw,yw],[x,y]]), 'euclidean')[0]
        return Q * k0(r/self.B)/(2*pi*self.K*self.b)                        # note: k0() = modified Bessel function of the first kind of zero order
    def Interface(self):                                                    # create a data input window for modifying hydro parameter values
        root = Tk()
        container = Frame(root)
        container.grid()
        Label(container,text='Hydraulic Properties',font=('Courier',14)).grid()        
        labels = ['Aquifer hydraulic conductivity (horizontal)','Aquifer thickness','Aquitard hydraulic conductivity (vertical)','Aquitard thickness','Aquifer porosity','Reference head']
        params = [self.K,self.b,self.K_c,self.b_c,self.phi,self.h0]
        entry = []
        for i in xrange(len(labels)):
            Label(container,text=labels[i]).grid(row=i+1,sticky=W)
            entry.append(Entry(container))
            entry[i].grid(row=i+1,column=1)
            entry[i].insert(0,str(params[i]))
        btn = Button(container,text='UPDATE',command=lambda:self.Button_click(entry)).grid(column=1)
        root.mainloop()
    def Button_click(self,entry):
        self.K = float(entry[0].get())
        self.b = float(entry[1].get())
        self.K_c = float(entry[2].get())        
        self.b_c = float(entry[3].get())
        self.phi = float(entry[4].get())        
        self.h0 = float(entry[5].get())        
        self.WriteParams()        
    def WriteParams(self):                                                
        # write present value set to text file
        output_file = open('hydro.txt','w')
        output_file.writelines(['K','\t',str(self.K),'\n'])
        output_file.writelines(['b','\t',str(self.b),'\n'])
        output_file.writelines(['K_c','\t',str(self.K_c),'\n'])
        output_file.writelines(['b_c','\t',str(self.b_c),'\n'])
        output_file.writelines(['phi','\t',str(self.phi),'\n'])        
        output_file.writelines(['h0','\t',str(self.h0),'\n'])        
        output_file.close()

class LinearSegment:
    def __init__(self,x0,y0,xf,yf):
        self.x0 = x0                                                        # starting and ending points
        self.y0 = y0
        self.xf = xf
        self.yf = yf
        self.x = (x0+xf)/2.                                                 # midpoint
        self.y = (y0+yf)/2.
        self.L = distance.pdist(array([[x0,y0],[xf,yf]]), 'euclidean')[0]   # length of segment
        self.xt0 = -0.5*self.L                                              # x-coordinates of segment endpoints in translated & rotated domain (integration endpoints)
        self.xtf = 0.5*self.L        
        if xf != x0:
            self.m = (self.yf-self.y0)/(self.xf-self.x0)                    # slope
            self.theta = arctan(self.m)                                     # orientation of element with respect to coordinate system (for grid rotation)
            self.alpha = arctan(-1./self.m)                                 # orientation of element normal with respect to coordinate system (for local hydraulic gradient projection)
        else:
            self.theta = pi/2.
            self.alpha = 0.
    def Transform(self,x,y):                                                ### translate and rotate coordinate (x,y) so it corresponds to an orientation with element parallel to x-axis, with origin at element midpoint
        x_t = x - self.x                                                    # translate to offset from element midpoint
        y_t = y - self.y
        x_prime = x_t*cos(-self.theta) - y_t*sin(-self.theta)
        y_prime = x_t*sin(-self.theta) + y_t*cos(-self.theta)
        return x_prime,y_prime

class FaultIntegrated(LinearSegment):                                       ##### container for collection of fault segments, used to spawn N Fault objects
    def __init__(self,x0,y0,xf,yf,R,N):
        LinearSegment.__init__(self,x0,y0,xf,yf)
        self.R = R                                                          # relative conductance (0 to 1)
        self.N = N                                                          # number of segments

class Fault(LinearSegment):                                                 ##### doublet element, integrated along a line segment
    def __init__(self,x0,y0,xf,yf,R):
        LinearSegment.__init__(self,x0,y0,xf,yf)
        self.R = R                                                          # relative conductance (0 to 1)
        self.R0 = 1.0                                                       # initial relative resistance (replaced by R after weight is set)
        self.S = 1.0                                                        # weighting factor; dummy initial value (used to define GW flow balance matrix)
        self.xc = self.x + cos(self.alpha)*dl                               # compliance point,normal to element where strength term is set
        self.yc = self.y + sin(self.alpha)*dl
    def DeltaH(self,x,y,hydro):                                             ### return head perturbation associated with doublet element, integrated along a finite length, at coordinate (x,y)
        x_prime,y_prime = self.Transform(x,y)
        return self.R0 * self.S * (arctan2(-y_prime,self.xt0-x_prime) - arctan2(-y_prime,self.xtf-x_prime))             # note that the integral of a point doublet cos(theta)/r is atan(x/y)
    def DirGrad(self,element,hydro):                                        ### return the projection of the hydraulic gradient vector along direction alpha associated with 'element'
        grad_x = (element.DeltaH(self.xc+dx,self.yc,hydro) - element.DeltaH(self.xc-dx,self.yc,hydro))/(2.*dx)
        grad_y = (element.DeltaH(self.xc,self.yc+dx,hydro) - element.DeltaH(self.xc,self.yc-dx,hydro))/(2.*dx)  
        return -dot([grad_x,grad_y],[cos(self.alpha),sin(self.alpha)])  

class HorizInjector(LinearSegment):                                         ##### finite-length injector, created by integrating a well element along a line segment
    def __init__(self,x0,y0,xf,yf,Qx):
        LinearSegment.__init__(self,x0,y0,xf,yf)        
        self.Qx = Qx                                                        # volumetric injection rate per unit length
    def DeltaH(self,x,y,hydro):                                             ### return head perturbation from injector by numerically integrating the line source along its length, with respect to coordinate (x,y)
        x_prime,y_prime = self.Transform(x,y)
        return quad(hydro.LeakyPointSource,self.xt0,self.xtf,args=(0.0,x_prime,y_prime,self.Qx))[0]

class Well:                                                                 ##### point source/sink
    def __init__(self,x0,y0,Q):
        self.x0 = x0
        self.y0 = y0
        self.Q = Q
    def DeltaH(self,x,y,hydro):                                             ### return head perturbation associated with well operation at coordinate (x,y)
        return hydro.LeakyPointSource(self.x0,self.y0,x,y,self.Q)


##########################################
#
# script
#
##########################################


def AEM():

    hydro = Hydro()
    print 'Read aquifer properties.'

    base = Base()
    print 'Read grid specifications and regional GW surface characteristics.'

    well,injector,fault_int = ReadElements()                    # read elements file to inform user interface (if modifications are to be done)
    interface = ElementInterface(base,well,injector,fault_int)
    well,injector,fault_int = ReadElements()                    # re-read (posisbly modified-) elements file to inform model
    print 'Read analytic elements.'

    # spawn fault segments from integrated fault(s) and assign strength terms
    for i in xrange(len(fault_int)):
        fault = []
        del_x = (fault_int[i].xf - fault_int[i].x0)/fault_int[i].N
        del_y = (fault_int[i].yf - fault_int[i].y0)/fault_int[i].N                
        for j in xrange(fault_int[i].N):
            fault.append(Fault(fault_int[i].x0+j*del_x,fault_int[i].y0+j*del_y,fault_int[i].x0+(j+1)*del_x,fault_int[i].y0+(j+1)*del_y,fault_int[i].R))
        print 'Spawned fault segments.'
    if len(fault_int):
        A,B = WeightMatrix(well,injector,fault,hydro)
        if len(fault) == 1: x = B/A[0]
        else: x=linalg.solve(A,B)
        for i in xrange(len(fault)):
            fault[i].S = x[i]               # assign weights
            fault[i].R0 = fault[i].R        # assign hydraulic resistances
        print 'Calculated strength terms for fault elements.'
    else:
        fault = []

    # particle tracking post-processor
    track = Track()
    if track.track:
        print 'Performing particle tracking ...'
        trail = track.MoveODE(base,well,injector,fault,hydro)

    print 'Populating grid points, plotting, and writing output ...'
    base.ProcessOutput(well,injector,fault,hydro,trail)

    print 'Done.'


############ run script

AEM()

