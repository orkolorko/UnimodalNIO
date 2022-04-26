alphas = 2:0.1:8
epsilons = 0:0.01:1.0

begin 
    n = 10000
    n_exp = 200

    function T_alpha(alpha,x)
        return 1-2*abs(x)^alpha
    end

    function T_prime_alpha(alpha,x)
        if x < 0
            return 2*alpha*abs(x)^(alpha-1)
        else
            return 2*alpha*(x)^(alpha-1)
        end
    end

    function periodic(x::Real)
        if x<-1
            return 1+(x+1)
        elseif x>1
            return -2+x
        else
            return x
        end
    end

    function compute_orbit_with_noise(n, epsilon, T)
        noise = epsilon*rand(n).-epsilon/2
        x_0 = rand(1)[1]
        orbit = zeros(n)
        orbit[1] = x_0+noise[1]
        for i in 2:n
            orbit[i]= periodic(T(orbit[i-1]) + noise[i])
        end
        return orbit
    end

    function compute_lyap(orbit, T_prime)
        sum = 0
        for i in 1:length(orbit)
            sum += log(abs(T_prime(orbit[i])))
        end
        return sum/length(orbit)
    end


    function compute_avg_lyap(epsilon, n, n_exp, T, T_prime)
        avg_lyap = 0
        for i in 1:n_exp
            orbit = compute_orbit_with_noise(n, epsilon, T)
            lyap = compute_lyap(orbit, T_prime)
            if lyap!= Inf
                avg_lyap+=lyap
            end
        end
        return avg_lyap/=n_exp
    end
end

lyapgrid = Array{Float64,2}(undef, length(alphas), length(epsilons))

for i in 1:length(alphas)
    println(alphas[i])
    T = x-> T_alpha(alphas[i], x)
    T_prime = x-> T_prime_alpha(alphas[i], x)
    for j in 1:length(epsilons)
        lyapgrid[i, j] = compute_avg_lyap(epsilons[j], n, n_exp, T, T_prime)
    end
end

using JLD
save("lyapgrid.jld", "lyapunov_grid", lyapgrid)
