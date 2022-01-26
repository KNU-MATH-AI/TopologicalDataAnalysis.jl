###################################정렬함수####################################
###1번###m_1이 작으면 M[:,n_1]을 제일왼쪽으로
###2번###m_2가 작으면 M[n_2,:]을 제일위로
###3번-1### 혹시 1열이 모두 첫원소의 배수이면 제일 마지막열에 이 수 곱하고 1열 나누기.그러면 1만들어짐
###3번-2### 3-1아니면 mod로 행과열 모두 나누기

###4번### 다시 1번부터
function eli_zoro_row_col(M::Matrix{T}) where T <: Integer 
    m,n = size(M)
    C_1 = []
    R_1 = []
    for i_1 in 1:m
        if count(M[i_1,:] .!= 0) == 0
            push!(R_1, i_1)
        end
    end      
    for i_2 in 1:n
        if count(M[:,i_2] .!= 0) == 0
            push!(C_1, i_2)
        end
    end      
    M = M[1:end .∉ [R_1], 1:end .∉ [C_1] ]  
    C_1 = []
    R_1 = []
    return M                       
end                                                     #zero행 또는 열을 없애준다
function cswap!(M::Matrix{T}, j1::T,j2::T) where T <: Integer
    j1 == j2 && return M
    M[:,[j1,j2]] = M[:,[j2,j1]]
    return M
end
function rswap!(M::Matrix{T}, i1::T,i2::T) where T <: Integer
    i1 == i2 && return M
    M[[i1,i2],:] = M[[i2,i1],:]
    return M
end
function zero_inf(x::Vector{T}) where T <: Integer
    if x > 0
        return x
    elseif x < 0
        return -x
    else
        return Inf
    end
end
function nonzero_min(M::Vector{T}) where T <: Integer
    nonzero_M = zero_inf.(M)
    m,n = findmin(nonzero_M)  ##최소숫자,위치찾기
    return M[n],n
end
function move_min_1_1(M::Matrix{T}) where T <: Integer  
    m_1,n_1 = nonzero_min(M[1,:])    ##첫행에서 (절대값이)가장 작은 수,위치
    m_2,n_2 = nonzero_min(M[:,1])    ##첫열에서 (절대값이)가장 작은 수,위치
    if abs(m_1) <= abs(m_2)
        cswap!(M, 1,n_1)       ##첫열의 (절대값이)최소수가 가장작으면 
    else
        rswap!(M, 1,n_2)       ##첫행의 (절대값이)최소수가 가장작으면 
    end    
end    
function isdivisible(M::Matrix{T}) where T <: Integer 
    m,n = size(M)
    for i in 2:m
        for j in 2:n
            if Int(M[i,j] % M[1,1]) != 0
                return false
            else 
                return true   
            end
        end
    end            
end    
function C_modular(M::Matrix{T}) where T <: Integer 
    a = M[1,2:end] .% M[1,1]                                  ##행렬에서 나누기
    s = size(a, 1)
    for i in 1:s
        if a[i] < 0
            if M[1,1] >0
                a[i] += M[1,1]
            else
                a[i] -= M[1,1] 
            end    
        end    
    end    
    b = convert(Vector{Int64},( M[1,2:end] .- a ) ./ M[1,1] )  ##행렬에서 몫 구하기
    M[:,2:end] = M[:,2:end] - (b)' .*M[:,1]                   ##column eli 
    return M
end    
function R_modular(M::Matrix{T}) where T <: Integer 
    a = M[2:end,1] .% M[1,1]                                  ##행렬에서 나누기
    s = size(a, 1)
    for i in 1:s
        if a[i] < 0
            if M[1,1] >0
                a[i] += M[1,1]
            else
                a[i] -= M[1,1] 
            end    
        end    
    end    
    b = convert(Vector{Int64},( M[2:end,1] .- a ) ./ M[1,1])  ##행렬에서 몫 구하기
    M[2:end,:] = M[2:end,:] - (b) .*M[1,:]'                   ##row eli 
    return M
end    
###########################################################################
function makeone(M::Matrix{T}) where T <: Integer 
    move_min_1_1(M) 
    C_modular(M)
    R_modular(M)       ##여기까지면 1,1을 제외하고 모두 0
    if count(M[1,2:1:end] .!= 0) == 0 && count(M[2:1:end,1] .!= 0) == 0
        if M[1,1] == 1
            push!(M_1, M[1,1])
            return M
        elseif M[1,1] == -1
            push!(M_1, -M[1,1])
            return M
        elseif isdivisible(M) == true
            push!( M_1, abs(M[1,1]) )
            return M                            ##이 숫자가 모두 나눠지면 return M      
        else 
            M[1,:]= M[1,:] .+M[2,:]             ##2행을 1행에 더하기
            return makeone(M)
        end
    end
    return makeone(M)
end
###########################################################################
function end_row(M::Matrix{T}) where T <: Integer 
    move_min_1_1(M) 
    C_modular(M)
    if count(M[1,2:1:end] .!= 0) == 0
        if M[1,1] > 0
            push!(M_1, M[1,1])
        elseif M[1,1] < 0
            push!(M_1, -M[1,1])  
        end      
    else  
        return end_row(M)
    end    
end    
#############################S_N_F계산함수#######################################
M_1 = []
function S_N_F(M::Matrix{T}) where T <: Integer 
    eli_zoro_row_col(M)            ##0행0열을 걸러낸다.
    makeone(M)
    M = M[2:end,2:end]             #M_1에 M의 1행1열 숫자 남기고 1행1열부분은 버린다.
    if size(M, 1) == 1
        end_row(M)
        return M_1
    end                             #1행만 남았으면 그만 
    return  S_N_F(M)
end






#########################
#########################
################################################################
################################################################
function eli_zoro_row_col(M::Matrix{T}) where T <: Integer 
    m,n = size(M)
    C_1 = []
    R_1 = []
    for i_1 in 1:m
        if count(M[i_1,:] .!= 0) == 0
            push!(R_1, i_1)
        end
    end      
    for i_2 in 1:n
        if count(M[:,i_2] .!= 0) == 0
            push!(C_1, i_2)
        end
    end      
    M = M[1:end .∉ [R_1], 1:end .∉ [C_1] ]  
    C_1 = []
    R_1 = []
    return M                       
end                                                     #zero행 또는 열을 없애준다
function cswap!(M::Matrix{T}, j1::T,j2::T) where T <: Integer
    j1 == j2 && return M
    M[:,[j1,j2]] = M[:,[j2,j1]]
    return M
end
function rswap!(M::Matrix{T}, i1::T,i2::T) where T <: Integer
    i1 == i2 && return M
    M[[i1,i2],:] = M[[i2,i1],:]
    return M
end
function zero_inf(x::Real)
    if x > 0
        return x
    elseif x < 0
        return -x
    else
        return Inf
    end
end
function nonzero_min(M::Vector{T}) where T <: Integer
    nonzero_M = zero_inf.(M)
    m,n = findmin(nonzero_M)  ##최소숫자,위치찾기
    return M[n],n
end
function move_min_1_1(M)      
    m_1,n_1 = nonzero_min(M[1,:])    ##첫행에서 (절대값이)가장 작은 수,위치
    m_2,n_2 = nonzero_min(M[:,1])    ##첫열에서 (절대값이)가장 작은 수,위치
    if abs(m_1) <= abs(m_2)
        cswap!(M, 1,n_1)       ##첫열의 (절대값이)최소수가 가장작으면 
    else
        rswap!(M, 1,n_2)       ##첫행의 (절대값이)최소수가 가장작으면 
    end    
end    
function isdivisible(M) 
    m,n = size(M)
    for i in 2:m
        for j in 2:n
            if Int(M[i,j] % M[1,1]) != 0
                return false
            else 
                return true   
            end
        end
    end            
end    
function C_modular(M)
    a = M[1,2:end] .% M[1,1]                                  ##행렬에서 나누기
    s = size(a, 1)
    for i in 1:s
        if a[i] < 0
            if M[1,1] >0
                a[i] += M[1,1]
            else
                a[i] -= M[1,1] 
            end    
        end    
    end    
    b = convert(Vector{Int64},( M[1,2:end] .- a ) ./ M[1,1] )  ##행렬에서 몫 구하기
    M[:,2:end] = M[:,2:end] - (b)' .*M[:,1]                   ##column eli 
    return M
end    
function R_modular(M)
    a = M[2:end,1] .% M[1,1]                                  ##행렬에서 나누기
    s = size(a, 1)
    for i in 1:s
        if a[i] < 0
            if M[1,1] >0
                a[i] += M[1,1]
            else
                a[i] -= M[1,1] 
            end    
        end    
    end    
    b = convert(Vector{Int64},( M[2:end,1] .- a ) ./ M[1,1])  ##행렬에서 몫 구하기
    M[2:end,:] = M[2:end,:] - (b) .*M[1,:]'                   ##row eli 
    return M
end    
###########################################################################
function makeone(M)
    move_min_1_1(M) 
    C_modular(M)
    R_modular(M)       ##여기까지면 1,1을 제외하고 모두 0
    if count(M[1,2:1:end] .!= 0) == 0 && count(M[2:1:end,1] .!= 0) == 0
        if M[1,1] == 1
            push!(M_1, M[1,1])
            return M
        elseif M[1,1] == -1
            push!(M_1, -M[1,1])
            return M
        elseif isdivisible(M) == true
            push!( M_1, abs(M[1,1]) )
            return M                            ##이 숫자가 모두 나눠지면 return M      
        else 
            M[1,:]= M[1,:] .+M[2,:]             ##2행을 1행에 더하기
            return makeone(M)
        end
    end
    return makeone(M)
end
###########################################################################
function end_row(M::Matrix{T}) where T <: Integer 
    move_min_1_1(M) 
    C_modular(M)
    if count(M[1,2:1:end] .!= 0) == 0
        if M[1,1] > 0
            push!(M_1, M[1,1])
        elseif M[1,1] < 0
            push!(M_1, -M[1,1])  
        end      
    else  
        return end_row(M)
    end    
end    
#############################S_N_F계산함수#######################################
M_1 = []
function S_N_F(M)
    eli_zoro_row_col(M)            ##0행0열을 걸러낸다.
    makeone(M)
    M = M[2:end,2:end]              #M_1에 M의 1행1열 숫자 남기고 1행1열부분은 버린다.
    if size(M, 1) == 1
        end_row(M)
        M_2 = M_1
        global M_1 = []
        return M_2
    end                             #1행만 남았으면 그만 
    return  S_N_F(M)
end


######################################################################################
######################################################################################
######################################################################################