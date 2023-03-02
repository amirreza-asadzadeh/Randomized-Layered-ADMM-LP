using LinearAlgebra, StatsBase

function ADMMdecoder_RL(DecoderMethod, γ, maxIteration, ϵ, μ, ρ, α, N_ET, ϕ)

    ###
    # This function applies non-penalized/l2-penalized ADMM-LP decoding using
    # the randomized layered (RL) scheduling.
    #
    # DecoderMethod: "LP" for non-penalized, "L2" for l2-penalized decoder
    # γ: the Log Likelihood Ratio (LLR) vector of size N
    # maxIteration: the maximum allowed number of iterations
    # ϵ: stopping precision, usually set to 10^(-5)
    # μ: the Lagrangian penalty parameter, usually set to 3
    # ρ: the over-relaxation parameter, generally 1 < ρ < 2, ρ = 1 for non-over-relaxed case, recommended to set to 1.9
    # α: the parameter for l2-penalty term, usually set to 0.8
    # N_ET: the number of iteration to check early termination
    #
    # Returns:
    # x_dec: a binary (0/1) array of size N, the decoded output
    # iter: the number of iterations it took for the decoder to decode
    #
    ###

    λ = Array{Array{Float64, 1}, 1}(undef, M)
    zReplica = Array{Array{Float64, 1}, 1}(undef, M)
    zOld = Array{Array{Float64, 1}, 1}(undef, M)
    Pjx = Array{Array{Float64, 1}, 1}(undef, M)

    for j = 1:M
        λ[j] = zeros(CheckDegree[j])
        zReplica[j] = ones(CheckDegree[j]).*0.5
        zOld[j] = ones(CheckDegree[j]).*0.5
        Pjx[j] = ones(CheckDegree[j]).*0.5
    end

    dis_vec = zeros(M)
    x_dec = ones(N).*0.5
    x_hard = similar(x_dec)
    temp = zeros(N)

    iter = 1
    feas_1 = 1e3; #holds for sum over j of abs(P_jx^k - z_j^k)
    feas_2 = 1e3; #holds for sum over j of abs(z_j^k - z_j^(k-1))
    feas_tol = ϵ^2 * M * maximum(CheckDegree)

    count_msg = 0

    while( count_msg < maxIteration * M && (feas_1 >= feas_tol || feas_2 >= feas_tol) )

        if count_msg < M # first iteration
            feas_1 = 0
            feas_2 = 0

            for j = 1:M

                # Update Variable's Information: L_{i}
                for i in F_list[j]
                    temp[i] = 0
                    for k in G_list[i]
                        IDX  = findall(x -> x ==i, F_list[k])[1]
                        temp[i] += zReplica[k][IDX] - (λ[k][IDX] / μ)
                    end
                    temp[i] -= γ[i] / μ
                end

                # Update VN to CN messages: L_{i ->  j}
                if occursin("LP", DecoderMethod)
                    for i in F_list[j]
                        x_dec[i] = min( max( ( temp[i] ) / (VariableDegree[i]), 0.0), 1.0)
                    end
                elseif occursin("L2", DecoderMethod)
                    for i in F_list[j]
                        x_dec[i] = min( max( ( temp[i] - (α / μ) ) / ( VariableDegree[i] - ( 2 * α / μ ) ), 0.0), 1.0)
                    end
                end

                # Update CN to VN messages: L_{j -> i}
                zOld[j] .= zReplica[j]
                for k = 1:CheckDegree[j]
                    IDX = F_list[j][k]
                    Pjx[j][k] = x_dec[IDX]
                end
                zReplica[j] = parpolyproj( ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .+ ( λ[j] ./ μ ) )
                dis_vec[j] = norm(zReplica[j] .- 0.5)
                difference1 = ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .- zReplica[j]
                λ[j] .+= μ .* ( difference1 )
                feas_1 += norm( difference1 )^2
                difference2 = zReplica[j] - zOld[j]
                feas_2 += norm( difference2 )^2
                count_msg += 1
            end

        else
            feas_1 = 0
            feas_2 = 0

            for k = 1:M

                weights = Weights(exp.(-ϕ.*dis_vec))
                j = sample(1:M, weights) # the sampled check index

                # Update Variable's Information: L_{i}
                for i in F_list[j]
                    temp[i] = 0
                    for k in G_list[i]
                        IDX  = findall(x -> x ==i, F_list[k])[1]
                        temp[i] += zReplica[k][IDX] - (λ[k][IDX] / μ)
                    end
                    temp[i] -= γ[i] / μ
                end

                # Update VN to CN messages: L_{i ->  j}
                if occursin("LP", DecoderMethod)
                    for i in F_list[j]
                        x_dec[i] = min( max( ( temp[i] ) / (VariableDegree[i]), 0.0), 1.0)
                    end
                elseif occursin("L2", DecoderMethod)
                    for i in F_list[j]
                        x_dec[i] = min( max( ( temp[i] - (α / μ) ) / ( VariableDegree[i] - ( 2 * α / μ ) ), 0.0), 1.0)
                    end
                end

                # Update CN to VN messages: L_{j -> i}
                zOld[j] .= zReplica[j]
                for k = 1:CheckDegree[j]
                    IDX = F_list[j][k]
                    Pjx[j][k] = x_dec[IDX]
                end
                zReplica[j] = parpolyproj( ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .+ ( λ[j] ./ μ ) )
                dis_vec[j] = norm(zReplica[j] .- 0.5)
                difference1 = ( ρ .* Pjx[j] ) .+ ( ( 1 - ρ ) .* zOld[j] ) .- zReplica[j]
                λ[j] .+= μ .* ( difference1 )
                feas_1 += norm( difference1 )^2
                difference2 = zReplica[j] - zOld[j]
                feas_2 += norm( difference2 )^2
                count_msg += 1
            end
        end

        # Early Termination
        if count_msg % (N_ET * M) == 0
            HardDecode(x_dec, x_hard)
            sum_syn = 0
            for j = 1:M
                sum_syn += sum( x_hard[ F_list[j] ] ) % 2
            end
            if sum_syn == 0
                x_dec = x_hard
                break
            end
        end
    end

    iter = count_msg / M
    HardDecode(x_dec, x_hard)
    x_dec = x_hard
    return x_dec, iter
end


function HardDecode(x, x_hard)

    ###
    # This function assigns binary 0/1 values to a vector of continous values. If a value is above 0.5, it is set to 1,
    # otherwise, it is set to 0. The function updates x_hard in-place to enhance performance.
    #
    # x: input vector with continous variables
    # x_hard: the output binary vector, it should be allocated with the same size as x
    #
    # Returns:
    # nothing
    #
    ###

    for i =1:N
        if x[i] <= 0.5
            x_hard[i] = 0;
        else
            x_hard[i] = 1;
        end
    end
    return nothing
end
