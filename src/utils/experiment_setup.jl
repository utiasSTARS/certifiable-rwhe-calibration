include("./rotation_noise.jl")

function random_test_instance(κ, vt, mean_rot, rot_std_dev, mean_trans, trans_std_dev, T_gt, n_motions, langevin=true)
    A = []
    B = []
    for _ in 1:n_motions
        Ra = random_rotation_sample_normal_magnitude(mean_rot, rot_std_dev)
        ta = random_direction()*mean_trans + randn(3)*trans_std_dev
        Ai = [Ra ta; 0 0 0 1]
        push!(A, Ai)
        Bi = inv(T_gt)*Ai*T_gt
        if langevin
            Bi[1:3, 1:3] = Bi[1:3, 1:3]*langevin_sample(κ)
        else
            Bi[1:3, 1:3] = Bi[1:3, 1:3]*exp(so3(randn(3)*sqrt(convert_langevins_concentration_to_lie_variance(κ))))
        end
        Bi[1:3, 4] = Bi[1:3, 4] + randn(3)*sqrt(vt)
        push!(B, Bi)
    end
    return A, B
end

function vec_of_matrices_to_tensor(A)
    n = length(A)
    A_out = zeros(n, size(A[1], 1), size(A[1], 2))
    for i in 1:n
        A_out[i, :, :] = A[i]
    end
    return A_out
end