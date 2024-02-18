include("./euler_solver.jl")

using .euler_ode_solve
using PyPlot

const m = 1.0                                                                   # just a constant which will serve as central value for initializing masses
const k = 1.0                                                                   # just a constant which will serve as central value for initializing stiffnesses
const num_particles = 100_000                                                   # number of particles
const h = 5e-3                                                                  # step size for euler
const num_steps = Int(1e4)                                                      # the number of steps the simulation is to run for
const m_arr = [rand((m-m*0.2):(0.001):(m+m*0.2)) for _ in 1:num_particles]      # initialize the masses
const k_arr = [rand((k-k*0.2):(0.001):(k+k*0.2)) for _ in 1:num_particles]      # initialize the stiffnesses
const x0 = 9.5                                                                  # just a constant which will serve as central value for initializing positions
const p0 = 1e-1                                                                 # just a constant which will serve as central value for initializing momenta                                                                         

function main()
    
    # initialize positions and momenta
    x_arr = [rand((x0-0.5):(0.001):(x0+0.5)) for _ in 1:num_particles]
    p_arr = [rand((p0-0.5):(0.001):(p0+0.5)) for _ in 1:num_particles]


    # taking g=1 :)
    function Fp(x, k, m)
        return -m*sin(x)
    end

    function Fx(p, m)
        return p/m
    end

    # driver code
    for _ in 1:num_steps

        # get the positions and the momenta for all the particles for the next timestep
        x_arr, p_arr = euler_ode_solve.euler_solve(h=h, num_of_particles=num_particles, Fx=Fx, Fp=Fp, mass_arr=m_arr, stiffness_arr=k_arr, x_arr=x_arr, p_arr=p_arr)

        # plot
        scatter(x_arr, p_arr, alpha=0.3, color="black", s=0.28)
        xlabel("x", fontsize=20)
        ylabel("p", fontsize=20)

        xlim(-10, 30)
        ylim(-4, 4)

        pause(0.001)
        cla()
    end
    axis("equal")
    show()

end

main()
