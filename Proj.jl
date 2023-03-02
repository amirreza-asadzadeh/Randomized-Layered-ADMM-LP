function probsimplexproj(v)

    ###
    # A function for projection onto probability simplex in the centered framework.
    # Function projects the point and assigns the output value to input for performance
    # enhancement.
    #
    # v: input point to be projected
    #
    # Returns:
    # nothing
    ###

    rho = sort(v, rev = true);
    d = length(v);
    u = Vector{Float64}(undef, d);
    sum = 0;
    for i = 1:d
        sum += rho[i];
        u[i] = (1/i) * (sum - 1);
    end
    i_star = 1;
    while( i_star <= d )
        if rho[i_star] > u[i_star]
            i_star += 1;
        else
            break
        end
    end
    i_star -= 1;
    for i = 1:d
        v[i] = max( v[i] - u[i_star] - 1/2, -1/2 );
    end
    return nothing
end


function parpolyproj(v)

    ###
    # This function projects a point onto the parity polytope. The input vector should
    # be given in the non-centered framework. Output is given in the non-centered framework
    # as well. However, the computations inside the function takes place in the centered
    # framework, as it is practically more efficient. The function updates the input vector
    # as the output vector to enhance performance of the code.
    #
    # v: input point to be projected
    #
    # Returns:
    # the projected point
    ###

    d = length(v);
    v -= ones(d) * (1/2); #shifting to fit in the centered projection framework
    f = zeros(d);
    for i = 1:d
        if v[i] >= 0
            f[i] = 1; # facet identifictation
        end
    end
    if sum(f)%2 == 0
        i_star = argmin( abs.(v) );
        f[i_star] = 1 - f[i_star];
    end
    v_tilde = similar(v);
    for i = 1:d
        v_tilde[i] = v[i] * (-1)^f[i]; # similarity transform
    end
    v_tildeclipped = similar(v);
    for i = 1:d
        v_tildeclipped[i] = max( min( v_tilde[i], 1/2 ), -1/2 ); # hypercube projeciton
    end
    if sum(v_tildeclipped) >= 1 - d/2 # membership test
      for i = 1:d
            v[i] = max( min( v[i], 1/2 ), -1/2 );
        end
    else
        probsimplexproj(v_tilde) # simplex projection
        for i = 1:d
            v[i] = v_tilde[i] * (-1)^f[i]; # similarity transform
        end
    end
    return v + ones(d) * (1/2) # invert shifting to fit in the non-centered projection framework
end
