# Compute and store statistics after each simulation
# Alex Tam, 18/07/2022

##### Data structure #####
"Data structure for statistical outputs of one simulation"
mutable struct Network_Statistics
    σ::Float64 # [pN/μm] Time-averaged mean normal stress
    Mean_Normal_Stress::Vector{Float64} # [pN/μm] Vector of mean normal stress (time-dependent)
    Stress::Vector{Vector{Float64}} # [pN/μm] Vector of stress components (time-dependent)
    Ia::Float64 # [-] Time-averaged aggregation statistic
    Ib::Float64 # [-] Time-averaged bundle statistic
    Ip::Float64 # [-] Time-averaged parallel bundle statistic
    Aggregation_Statistic::Vector{Float64} # [-] Aggregation index (time-dependent)
    Bundle_Statistic::Vector{Float64} # [-] Bundle formation index (time-dependent)
    Parallel_Statistic::Vector{Float64} # [-] Parallel bundle index (time-dependent)
    Bundle_Direction::Vector{Float64} # [rad] Principal direction of filaments in bundled structure
    # Time-dependent and time-averaged mean relative motor position
    mp::Float64 # [-] Time-averaged mean relative motor position
    Motor_Position::Vector{Float64} # [-] Mean relative motor position (time-dependent)
    # Time-dependent and time-averaged mean filament curvature
    K::Float64 # [-] Time-averaged mean integrated filament curvature
    Curvature::Vector{Float64} # [-] Mean relative motor position (time-dependent)
    # Time-dependent and time-averaged relative angle
    θ::Float64 # [-] Time-averaged mean motor angle
    Motor_Angle::Vector{Float64} # [-] Mean motor angle (time-dependent)
    # Time-dependent and time-averaged I2
    I2::Float64 # [-] Time-averaged two-filament index
    TwoFilamentIndex::Vector{Float64} # [-] Two-filament index (time-dependent)
end

##### Store network (time-averaged) statistics #####
"Store network statistics"
function network_statistics(parN, Force, Aggregation_Statistic, Bundle_Statistic, Parallel_Statistic, Bundle_Direction, Motor_Position, Curvature, Motor_Angle, TwoFilamentIndex, Lxx, Lyy)
    NetStats = Network_Statistics(0.0, Vector{Vector}(), Vector{Vector{Vector}}(), 0.0, 0.0, 0.0, Aggregation_Statistic, Bundle_Statistic, Parallel_Statistic, Bundle_Direction, 0.0, Motor_Position, 0.0, Curvature, 0.0, Motor_Angle, 0.0, TwoFilamentIndex) # Pre-allocate output structure
    times = (0:parN.nT-1).*parN.dt # Vector of times at which we obtain measurements
    NetStats.Ia = trapz(Aggregation_Statistic, times)/(times[end]) # Time-averaged aggregation statistic
    NetStats.Ib = trapz(Bundle_Statistic, times)/(times[end]) # Time-averaged bundle statistic
    NetStats.Ip = trapz(Parallel_Statistic, times)/(times[end]) # Time-averaged bundle statistic
    NetStats.mp = trapz(Motor_Position, times)/(times[end]) # Time-averaged mean relative motor position
    NetStats.K = trapz(Curvature, times)/(times[end]) # Time-averaged mean integrated filament curvature
    NetStats.θ = trapz(Motor_Angle, times)/(times[end]) # Time-averaged mean integrated filament curvature
    NetStats.I2 = trapz(TwoFilamentIndex, times)/(times[end]) # Time-averaged mean integrated filament curvature
    NetStats.Stress, NetStats.Mean_Normal_Stress, NetStats.σ = stress(parN, Force, Lxx, Lyy) # Obtain stress results
    return NetStats
end

##### Functions for general operations #####
"Compute time-integrated quantities using the trapezoidal rule"
function trapz(y, x)
    int = 0.0 # Pre-allocate integral
    for i = 1:length(y)-1
        int += 0.5*(y[i] + y[i+1])*(x[i+1]-x[i]); # Composite trapezoidal rule
    end
    return int
end

"Calculate a moving average"
function movingaverage(x::Vector, numofele::Int)
    BackDelta = div(numofele,2)
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    len = length(x)
    y = similar(x)
    for n = 1:len
        lo = max(1, n - BackDelta)
        hi = min(len, n + ForwardDelta)
        y[n] = mean(x[lo:hi])
    end
    return y
end

"Compute Euclidean distance between two nodes"
function euclid_distance(n1x::Float64, n1y::Float64, n2x::Float64, n2y::Float64)
    return sqrt((n2x - n1x)^2 + (n2y - n1y)^2)
end

##### Compute specific statistics #####
"Compute stress components from a vector of force vectors"
function stress(parN::Numerical_Parameters, Force::Vector{Vector{Float64}}, Lxx::Float64, Lyy::Float64)
    times = (0:parN.nT-1).*parN.dt # Vector of times at which we obtain measurements
    Stress = [[0.0, 0.0, 0.0, 0.0] for idx in 1:length(Force)] # [pN/um] Vector of stress components (σxx, σxy, σyx, σyy) for all time steps
    Mean_Normal_Stress = [0.0 for idx in 1:length(Force)] # [pN/um] Vector of mean normal stress for all time steps
    for i in eachindex(Force) # Loop over time steps
        Stress[i][1] = Force[i][1]/Lyy # σxx
        Stress[i][2] = Force[i][2]/Lyy # σxy
        Stress[i][3] = Force[i][3]/Lxx # σyx
        Stress[i][4] = Force[i][4]/Lxx # σyy
        Mean_Normal_Stress[i] = 0.5*(Stress[i][1] + Stress[i][4]) # Mean normal stress σ = 0.5*tr(Stress)
    end
    σ = trapz(Mean_Normal_Stress, times)/(times[end]) # Time-averaged mean normal stress
    return Stress, Mean_Normal_Stress, σ
end

"Determine whether two filaments are aggregated, bundled, parallel-bundled"
function bundled(s, f1::Actin_Filament, f2::Actin_Filament, dc::Float64, θc::Float64)
    # f1 and f2 are two filaments, dist is the distance threshold, θ is the angle threshold (radians)
    aggregate::Bool = false # Default outcome (no aggregation)
    bundle::Bool = false # Default outcome (no bundle)
    parallel::Bool = false # Default outcome (no parallel bundle)
    angles = Vector{Float64}() # Preallocate empty vector of angles
    for seg1 in f1.segments # Loop over segments
        for seg2 in f2.segments # Loop over segments
            for i in eachindex(seg1.t) # Loop over translations for periodic BCs
                for j in eachindex(seg2.t) # Loop over translations for periodic BCs
                    m1x, m1y, p1x, p1y = get_segment_nodes(f1, seg1, s, seg1.t[i]) # Get dimenionless, translated node positions
                    m2x, m2y, p2x, p2y = get_segment_nodes(f2, seg2, s, seg2.t[j]) # Get dimenionless, translated node positions
                    # 1. Check proximity of both minus ends
                    if euclid_distance(m1x, m1y, m2x, m2y) < dc
                        aggregate = true
                        # If sufficiently close, check relative orientation
                        θs1 = atan(p1y-m1y, p1x-m1x) # Segment 1 orientation, [-π,π]
                        θs2 = atan(p2y-m2y, p2x-m2x) # Segment 2 orientation, [-π,π]
                        if (abs(θs1-θs2) < θc) || (abs(θs1-θs2) > (2*π - θc)) # Orientations close (parallel)
                            bundle = true
                            parallel = true
                            push!(angles, θs1)
                            push!(angles, θs2)
                        end
                        if (abs(θs1-θs2) > (π-θc)) && (abs(θs1-θs2) < (π+θc)) # Orientations close (anti-parallel)
                            bundle = true
                            push!(angles, θs1)
                            push!(angles, θs2)
                        end
                    end
                    # 2. Check proximity of first minus end and second plus end
                    if euclid_distance(m1x, m1y, p2x, p2y) < dc
                        aggregate = true
                        # If sufficiently close, check relative orientation
                        θs1 = atan(p1y-m1y, p1x-m1x) # Segment 1 orientation, [-π,π]
                        θs2 = atan(p2y-m2y, p2x-m2x) # Segment 2 orientation, [-π,π]
                        if (abs(θs1-θs2) < θc) || (abs(θs1-θs2) > (2*π - θc)) # Orientations close (parallel)
                            bundle = true
                            parallel = true
                            push!(angles, θs1)
                            push!(angles, θs2)
                        end
                        if (abs(θs1-θs2) > (π-θc)) && (abs(θs1-θs2) < (π+θc)) # Orientations close (anti-parallel)
                            bundle = true
                            push!(angles, θs1)
                            push!(angles, θs2)
                        end
                    end
                    # 3. Check proximity of first plus end and second minus end
                    if euclid_distance(p1x, p1y, m2x, m2y) < dc
                        aggregate = true
                        # If sufficiently close, check relative orientation
                        θs1 = atan(p1y-m1y, p1x-m1x) # Segment 1 orientation, [-π,π]
                        θs2 = atan(p2y-m2y, p2x-m2x) # Segment 2 orientation, [-π,π]
                        if (abs(θs1-θs2) < θc) || (abs(θs1-θs2) > (2*π - θc)) # Orientations close (parallel)
                            bundle = true
                            parallel = true
                            push!(angles, θs1)
                            push!(angles, θs2)
                        end
                        if (abs(θs1-θs2) > (π-θc)) && (abs(θs1-θs2) < (π+θc)) # Orientations close (anti-parallel)
                            bundle = true
                            push!(angles, θs1)
                            push!(angles, θs2)
                        end
                    end
                    # 4. Check proximity of both plus ends
                    if euclid_distance(p1x, p1y, p2x, p2y) < dc
                        aggregate = true
                        # If sufficiently close, check relative orientation
                        θs1 = atan(p1y-m1y, p1x-m1x) # Segment 1 orientation, [-π,π]
                        θs2 = atan(p2y-m2y, p2x-m2x) # Segment 2 orientation, [-π,π]
                        if (abs(θs1-θs2) < θc) || (abs(θs1-θs2) > (2*π - θc)) # Orientations close (parallel)
                            bundle = true
                            parallel = true
                            push!(angles, θs1)
                            push!(angles, θs2)
                        end
                        if (abs(θs1-θs2) > (π-θc)) && (abs(θs1-θs2) < (π+θc)) # Orientations close (anti-parallel)
                            bundle = true
                            push!(angles, θs1)
                            push!(angles, θs2)
                        end
                    end
                end
            end
        end
    end
    return aggregate, bundle, parallel, angles
end

"Compute bundle statistic and bundle directions at one instant"
function bundle_statistic(state, af)
    # Pre-allocate
    total_count::Int = 0 # Initialising total number of filament pairs 
    aggregated_count::Int = 0 # Initialise number of aggregated filaments pairs
    bundled_count::Int = 0 # Initialise number of filament pairs in bundle
    parallel_count::Int = 0 # Initialise number of filaments in parallel bundle
    bundled_directions = Vector{Float64}() # Pre-allocate vector of bundled directions

    # Threshold parameters
    dc = 0.1
    θc = π/9

    # Loop over all filament pairs
    for i in eachindex(af)
        for j in eachindex(af)
            if j > i # Prevent double counting
                total_count += 1 # Total pair count
                aggregate, bundle, parallel, angles = bundled(state, af[i], af[j], dc::Float64, θc::Float64)
                if aggregate == true
                    aggregated_count += 1 # Aggregate pair count
                end
                if bundle == true
                    bundled_count += 1 # Bundled pair count
                    for a in angles
                        push!(bundled_directions, a)
                    end
                end
                if parallel == true
                    parallel_count += 1 # Parallel bundle count
                end
            end
        end
    end
    return aggregated_count/total_count, bundled_count/total_count, parallel_count/total_count, bundled_directions # Ia, Ib, Ip, all directions
end

"Mean (integrated) filament curvature at one instant"
function curvature(af, s, Lxx, Lxy, Lyx, Lyy)
    filament_curvatures = Vector{Float64}() # Pre-allocate mean curvature data
    # Loop over filaments
    for f in af
        seg_lengths = get_segment_lengths(f, s, Lxx, Lxy, Lyx, Lyy) # Dimensional segment lengths
        total_filament_curvature = 0.0 # Pre-allocate integrated curvature of current filament
        # Compute curvature at interior nodes
        if length(s.an[f.index]) >= 3
            for i = 1:length(s.an[f.index])-2 
                # Extract physical positions of plus node
                xp = s.an[f.index][i+2][1]*Lxx + s.an[f.index][i+2][2]*Lyx
                yp = s.an[f.index][i+2][1]*Lxy + s.an[f.index][i+2][2]*Lyy
                # Extract physical positions of centre node
                xc = s.an[f.index][i+1][1]*Lxx + s.an[f.index][i+1][2]*Lyx
                yc = s.an[f.index][i+1][1]*Lxy + s.an[f.index][i+1][2]*Lyy
                # Extract physical positions of minus node
                xm = s.an[f.index][i][1]*Lxx + s.an[f.index][i][2]*Lyx
                ym = s.an[f.index][i][1]*Lxy + s.an[f.index][i][2]*Lyy
                # Compute numerical second derivatives
                L01 = seg_lengths[i]
                L12 = seg_lengths[i+1] # Extract segment lengths
                avl = (L01 + L12)/2 # Average length
                ddfx = ((xp-xc)/L12 - (xc-xm)/L01)/avl # Numerical second derivative (x)
                ddfy = ((yp-yc)/L12 - (yc-ym)/L01)/avl # Numerical second derivative (y)
                # Add contribution to integrated curvature
                total_filament_curvature += sqrt(ddfx^2 + ddfy^2)*avl
            end
        end
        push!(filament_curvatures, total_filament_curvature)
    end
    return mean(filament_curvatures)
end

"Mean relative motor position at one instant"
function mean_relative_mp(s::State)
    # Pre-allocate
    motor_pos = Vector{Float64}()
    for i = 1:length(s.mp) # Loop over motors
        push!(motor_pos, s.mp[i][1]) # Store motor position on 1st filament
        push!(motor_pos, s.mp[i][2]) # Store motor position on 2nd filament
    end
    return mean(motor_pos)
end

"Mean angle between two filaments and two-filament index at an instant"
function two_filament_index(mm, s, parA, Lxx, Lxy, Lyx, Lyy)
    motor_angles = Vector{Float64}() # Pre-allocate vector of motor angles
    index = Vector{Float64}() # Pre-allocate index data
    # Loop over motors
    for i in eachindex(mm)
        m::Myosin_Motor = mm[i] # Extract current motor
        f1::Actin_Filament = m.f1
        f2::Actin_Filament = m.f2 # Extract filaments attached to current motor
        # Calculate filament lengths
        L1 = sum(get_segment_lengths(f1, s, Lxx, Lxy, Lyx, Lyy))
        L2 = sum(get_segment_lengths(m.f2, s, Lxx, Lxy, Lyx, Lyy))
        # Calculate motor positions
        x1, y1, x2, y2 = get_motor_pos(m, s, parA, Lxx, Lxy, Lyx, Lyy)
        m1x = (x1 - m.t1[1])*Lxx + (y1 - m.t1[2])*Lyx
        m1y = (x1 - m.t1[1])*Lxy + (y1 - m.t1[2])*Lyy
        m2x = (x2 - m.t2[1])*Lxx + (y2 - m.t2[2])*Lyx
        m2y = (x2 - m.t2[1])*Lxy + (y2 - m.t2[2])*Lyy
        # Find segments on which motor attaches
        seg1, seg2, r1, r2 = get_motor_relative_pos_segment(m, s, parA, Lxx, Lxy, Lyx, Lyy)
        # Calculate segment positions
        mx1, my1, px1, py1 = get_segment_nodes(f1, f1.segments[seg1], s, m.t1, Lxx, Lxy, Lyx, Lyy)
        mx2, my2, px2, py2 = get_segment_nodes(f2, f2.segments[seg2], s, m.t2, Lxx, Lxy, Lyx, Lyy)
        # Calculate angle between vectors
        vec1 = [px1 - mx1, py1 - my1] # Vector between motor and plus end of filament 1
        vec2 = [px2 - mx2, py2 - my2] # Vector between motor and plus end of filament 2
        cosine = (vec1[1]*vec2[1] + vec1[2]*vec2[2])/( sqrt(vec1[1]^2 + vec1[2]^2)*sqrt(vec2[1]^2 + vec2[2]^2)) # Cosine of angle between vectors
        # Prevent spurious errors that occur in rare instances
        if abs(cosine) <= 1
            theta = acos(cosine) # Angle between vectors
        elseif cosine > 1
            theta = acos(1) # Angle between vectors
        else
            theta = acos(-1) # Angle between vectors
        end
        push!(motor_angles, theta)
        # Compute two-filament index
        if all(s.mp[i] .<= 1)
            Id = ( 2*(L1*s.mp[i][1] + L2*s.mp[i][2])/(L1+L2) - 1)*sin(theta/2)^2 # Heuristic index based on four filament branches
        else
            Id = 0 # Return zero if motor is detached
        end
        push!(index, Id)
    end
    return mean(motor_angles), mean(index)
end