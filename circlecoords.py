import numpy as np

def circle(dim, plane, radius, n_points, *origin):
    theta = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    ones = np.ones_like(theta)
    if dim == 3:
        x0,y0,z0 = origin
        if plane == 'z':
            a = x0 + radius*np.cos(theta)
            b = y0 + radius*np.sin(theta)
            c = z0*ones

        elif plane == 'y':
            a = x0 + radius*np.cos(theta)
            b = y0*ones
            c = z0 + radius*np.sin(theta)

        else:
            a = x0*ones
            b = y0 + radius*np.cos(theta)
            c = z0 + radius*np.sin(theta)
        return a,b,c

    else:
        x0,y0 = origin
        a = x0 + radius*np.cos(theta)
        b = y0 + radius*np.sin(theta)
        return a,b