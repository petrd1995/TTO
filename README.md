    #Truss topology optimization
    Programs in this repo contain functions used to create and solve problems in truss topology optimization.
    
    The core program is TTO.py, where the class Truss is defined along with its various methods.


    #Initialization of class instance with variables defining the particular problem.
        The variables are as follows:
        
        **name** - name of the problem (for saving purposes)
        **x0** - length of the design domain in x direction
        **y0** - length of the design domain in y direction
        **z0** - length of the design domain in z direction
        **nx** - number of nodes in x direction
        **ny** - number of nodes in y direction
        **nz** - number of nodes in z direction
        **bc** - vector(np.array) containing boundary conditions, which are explained by example:
            
            bc = [[0,1,1],[3,1,0]]

            [node_nr, bcx, bcy, (bcz)]

            where two vectors mean we have two bcs (or bc in two nodes)
            the nodes in which the bcs are present are defined by the firs column of each
            vector, so the first bc will be set in node number 0 and the second in node
            number 3. The remaining numbers define if corresponding degree of freedom
            is constrained or not - 1 means constrained, 0 means free. Node 0 is therefore
            constrained in both x and y, whereas node 3 is only constrained in x, and 
            free in y. For 3D problems it is needed to add one more column for z direction.
        **F** - vector(np.array) containing forces, similar to bcs:

            F = [[2, 0, -1000, 0],[4, 10, -1000, 20]]
            
            [node_nr, Fx, Fy, (Fz)]

            where two vectors mean two forces are defined. First number in each vector (row)
            again means node number in which the force is present, and the other numbers
            mean components of the individual forces.
        **E** - Young's modulus E of the material (MPa)
        **r0** - initial radius of all bars (mm)
        **Vol0** - particular volume from which optimal design is to be found in further defined
            ratio, if Vol0 is to be defined by the volume of the design domain leave
            Vol0 = 0
        **ratio** - ratio of volume of final design to initial Vol0, that is ratio = Vol/Vol0
        **R0** - density of the material (kg/mm^3) for dynamic problems
        **kon** - konvergence criteria, which is later defined as difference of two consecutive
            design of areas of individual bars. kon = 1 means that when only improvement
            of one percent is achieved, the optimization loop ends

    *example of problem creation:*
        benchmark = Truss('benchmark', 100, 100, 100, 2, 2, 2, np.array([[0, 1, 1, 1], [1, 1, 1, 1],[2, 1, 1, 1], [3, 1, 1, 1]]), np.array([[4, 0, -1000, -600]]), 2.1e5, 5, 0, 0.2, 100, 1)

        to run the problem we have to call some methods from the class:
            benchmark.create_grid()
            benchmark.create_bars()
            benchmark.vec_len()
            benchmark.opt()

        
        