include("./euler_solver.jl")

using .euler_ode_solve
using PyPlot

const m = 0.4                                                               # just a constant which will serve as central value for initializing masses
const k = 1.0                                                               # just a constant which will serve as central value for initializing stiffnesses
const num_particles = 10_000                                                # number of particles
const h = 1e-2                                                              # step size of euler
const num_steps = Int(1e4)                                                  # number of timesteps for which the simulation has to run
const m_arr = [rand((m-0.02):(0.001):(m+0.02)) for _ in 1:num_particles]    # array for masses of all particles 
const k_arr = [rand((k-0.02):(0.001):(k+0.02)) for _ in 1:num_particles]    # array for stifness of oscillators
const x0 = 1.0                                                              # just a constant which will serve as central value for initializing positions
const p0 = 3.0                                                              # just a constant which will serve as central value for initializing momenta

function main()
    
    # initialise positions and momenta for all particles
    x_arr = [rand((x0-0.5):(0.001):(x0+0.5)) for _ in 1:num_particles]
    p_arr = [rand((p0-0.5):(0.001):(p0+0.5)) for _ in 1:num_particles]

    # m*dp/dt = Fp
    function Fp(x, k, m)
        return -0.5*k*x
    end

    # dx/dt = Fx
    function Fx(p, m)
        return p/m
    end

    # driver code
    for _ in 1:num_steps

        # solve the diff eqns of motion for the next timestep for all the particles
        x_arr, p_arr = euler_ode_solve.euler_solve(h=h, num_of_particles=num_particles, Fx=Fx, Fp=Fp, mass_arr=m_arr, stiffness_arr=k_arr, x_arr=x_arr, p_arr=p_arr)

        # plot
        scatter(x_arr, p_arr, alpha=1, color="black", s=0.01)
        xlabel("x", fontsize=20)
        ylabel("p", fontsize=20)

        xlim(-8, 8)
        ylim(-6, 6)

        pause(0.0001)
        cla()
    end
    

end

main()

