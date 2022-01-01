import numpy as np
 
def section_angle(chord_length, arclength, tol, max_iter,
                     min_theta=0,
                     max_theta=np.pi,
                     iters=0
                    ):
    """
    Recursive bisection search for the angle (and radius)
    of the circular section
    """
    theta = .5 * (max_theta + min_theta) # seems like a good guesstimate
    r = arclength / theta
    if iters > max_iter:
        print("max_iter reached in search for section angle")
        return theta, r
    halfchord = r * np.sin( theta / 2 )
    #print("chord estimate" , 2 * halfchord)
    #print("actual chord" , chord_length)
    #print("min_theta" , min_theta)
    #print("max_theta" , max_theta, '\n')
    if np.abs( halfchord - chord_length/2 ) < tol:
        return theta, r
    elif halfchord < chord_length/2:
        max_theta = theta
    else:
        min_theta = theta
    return section_angle(chord_length, arclength, 
                         min_theta=min_theta,
                         max_theta=max_theta,
                         iters=iters+1,
                         tol=tol,
                         max_iter=max_iter
                        )

def circular_section(chord_length, arclength, npts, tol=1e-3, max_iter=100):
    theta, r = section_angle(chord_length, arclength,
                             tol=tol,
                             max_iter=max_iter)
    angles = np.linspace(np.pi/2 - theta/2, np.pi/2 + theta/2, npts)
    complex_arc = r * np.exp( 1j * angles)
    x = np.real(complex_arc)
    y = np.imag(complex_arc)
    return x, y

"""
import matplotlib.pyplot as plt

def test():
    chord_length = 5
    arclength = 5.5
    npts = 100
    x, y = circular_section(chord_length, arclength, npts)
    y = y - y[0]
    x = x - x[0]
    z = x + 1j * y
    z = z * np.exp( -1j * np.pi/6 )
    plt.axis('equal')
    plt.plot(np.real(z), np.imag(z))
    print(y[0] - y[-1])
    
test()
"""