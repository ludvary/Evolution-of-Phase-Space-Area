module euler_ode_solve

    # making a function which takes in coords of many particles and spits out coords of those particles at next time instant

    function euler_solve(;h :: AbstractFloat, num_of_particles :: Integer, Fx :: Function, Fp :: Function, mass_arr :: Vector{Float64}, stiffness_arr :: Vector{Float64}, x_arr :: Vector{Float64}, p_arr :: Vector{Float64})

            for i in 1:num_of_particles
                x_arr[i] += h*Fx(p_arr[i], mass_arr[i]) 
                p_arr[i] += h*Fp(x_arr[i], stiffness_arr[i], mass_arr[i])
            end

           return x_arr, p_arr 
            
        end


    # making a function which takes in coord of one particle and returns coords for many next time instances

    function euler_solve_single_particle(;h :: AbstractFloat, num_timesteps :: Integer, Fx :: Function, Fp :: Function, mass :: AbstractFloat, stiffness :: AbstractFloat, x0 :: AbstractFloat, p0 :: AbstractFloat)

        x_arr = zeros(AbstractFloat, num_timesteps)
        p_arr = zeros(AbstractFloat, num_timesteps)
        x = x0
        p = p0

        for i in 1:num_timesteps
            x_arr[i] = x
            p_arr[i] = p
            x += h*Fx(p, mass)
            p += h*Fp(x, stiffness)
        end

        return x_arr, p_arr
    end



end
