############# R_1 & R_2 ###################################
function R_1(x,y,m,s,t)
    z = ((1-s)x+s*y)%m
    if z < 1
        z += m
    end    
    return z
end  
function R_2(x,y,m,s,t)
    z = ((1-t)y+t*x)%m
    if z < 1
        z += m
    end    
    return z
end  
#############################################################
#############5dim Bounmdary map##############################  
count_1 = 1
function bound_matrix_5(m,s,t)
    n_boundary = zeros(Int64, m^(4), m^(5))
    for i = 1:m
        for j = 1:m
            for k = 1:m
                for l = 1:m
                    for n = 1:m
                        abc_ = []
                        abc_ = (j, k, l, n)                                                                  ## 1 ## 
                        abc_ = [abc_; (-R_1(i,j,m,s,t),-R_1(i,k,m,s,t),-R_1(i,l,m,s,t),-R_1(i,n,m,s,t))]     ## 2 ##   
                        abc_ = [abc_; (-R_2(i,j,m,s,t), -k, -l, -n)]                                         ## 3 ##   
                        abc_ = [abc_; (i, R_1(j,k,m,s,t), R_1(j,l,m,s,t),R_1(j,n,m,s,t))]                    ## 4 ## 
                        abc_ = [abc_; (R_2(i,k,m,s,t), R_2(j,k,m,s,t), l, n)]                                ## 5 ## 
                        abc_ = [abc_; (-i, -j,-R_1(k,l,m,s,t),-R_1(k,n,m,s,t))]                              ## 6 ##   
                        abc_ = [abc_; (-R_2(i,l,m,s,t), -R_2(j,l,m,s,t), -R_2(k,l,m,s,t),-n)]                ## 7 ##   
                        abc_ = [abc_; (i, j, k, R_1(l,n,m,s,t))]                                             ## 8 ##   
                        abc_ = [abc_; (R_2(i,n,m,s,t), R_2(j,n,m,s,t), R_2(k,n,m,s,t), R_2(l,n,m,s,t))]      ## 9 ##   
                        abc_ = [abc_; (-i, -j, -k, -l)]                                                      ## 10 ##   
                        for i_1 = 1:m
                            for j_1 = 1:m
                                for k_1 = 1:m
                                    for l_1 = 1:m
                                        for n_1 = 1:10
                                            if (i_1, j_1, k_1, l_1) == abc_[n_1]
                                                n_boundary[( (m^3)*(i_1-1) ) + (m^2)*(j_1-1) + (m*(k_1-1)) + (l_1), count_1] += 1
                                            elseif ((-1)*(i_1), (-1)*(j_1),(-1)*(k_1),l_1) == abc_[n_1]   
                                                n_boundary[( (m^3)*(i_1-1) ) + (m^2)*(j_1-1) + (m*(k_1-1)) + (l_1), count_1] += -1 
                                            end
                                        end    
                                    end 
                                end     
                            end
                        end                          ##여기서 tuple의 원소를 보고 n_boundary에 update한다. i,j,k가 뭔지는 중요하지않다.세로열만 중요.
                        global count_1 += 1      ##count_1에 +1 해준다.##   
                    end
                end
            end         
        end
    end 
    global count_1 = 1
    return n_boundary                      
end
##############################################################
#############4dim Bounmdary map###############################  
count_1 = 1
function bound_matrix_4(m,s,t)
    n_boundary = zeros(Int64, m^(3), m^(4))
    for i = 1:m
        for j = 1:m
            for k = 1:m
                for l = 1:m
                    abc_ = []
                    abc_ = (j, k, l)                                                        ## 1 ## 
                    abc_ = [abc_; (-R_1(i,j,m,s,t),-R_1(i,k,m,s,t),-R_1(i,l,m,s,t))]     ## 2 ##   
                    abc_ = [abc_; (-R_2(i,j,m,s,t), -k, -l)]                              ## 3 ##   
                    abc_ = [abc_; (i, R_1(j,k,m,s,t), R_1(j,l,m,s,t))]                    ## 4 ## 
                    abc_ = [abc_; (R_2(i,k,m,s,t), R_2(j,k,m,s,t), l)]                    ## 5 ## 
                    abc_ = [abc_; (-i, -j,-R_1(k,l,m,s,t))]                               ## 6 ##   
                    abc_ = [abc_; (-R_2(i,l,m,s,t), -R_2(j,l,m,s,t), -R_2(k,l,m,s,t))]    ## 7 ##   
                    abc_ = [abc_; (i, j, k)]                                              ## 8 ##   

                    for i_1 = 1:m
                        for j_1 = 1:m
                            for k_1 = 1:m
                                for l_1 = 1:8
                                if (i_1, j_1,k_1) == abc_[l_1]
                                    n_boundary[( (m^2)*(i_1-1) ) + (m*(j_1-1)) + (k_1), count_1] += 1
                                elseif ((-1)*(i_1), (-1)*(j_1),(-1)*(k_1)) == abc_[l_1]   
                                    n_boundary[( (m^2)*(i_1-1) ) + (m*(j_1-1)) + (k_1), count_1] += -1 
                                end    
                                end 
                            end     
                        end
                    end                          ##여기서 tuple의 원소를 보고 n_boundary에 update한다. i,j,k가 뭔지는 중요하지않다.세로열만 중요.
                    global count_1 += 1      ##count_1에 +1 해준다.##   
                end
            end         
        end
    end 
    global count_1 = 1
    return n_boundary                      
end
###############################################################
#############3dim Bounmdary map################################
count_1 = 1
function bound_matrix_3(m,s,t)
    n_boundary = zeros(Int64, m^(2), m^3)
    for i = 1:m
        for j = 1:m
            for k = 1:m
                abc_ = []
                abc_ = (j,k)                                           ## 1 ## 
                abc_ = [abc_; (-R_1(i,j,m,s,t), -R_1(i,k,m,s,t)) ]     ## 2 ##   
                abc_ = [abc_; (-R_2(i,j,m,s,t), -k)]                   ## 3 ##   
                abc_ = [abc_; (i, R_1(j,k,m,s,t))]                     ## 4 ## 
                abc_ = [abc_; (R_2(i,k,m,s,t), R_2(j,k,m,s,t))]        ## 5 ## 
                abc_ = [abc_; (-i, -j)]                                ## 6 ##   

                for i_1 = 1:m
                    for j_1 = 1:m
                        for k_1 = 1:6
                        if (i_1, j_1) == abc_[k_1]
                            n_boundary[(m*(i_1-1))+(j_1), count_1] += 1
                        elseif ((-1)*(i_1), (-1)*(j_1)) == abc_[k_1]   
                            n_boundary[(m*(i_1-1))+(j_1), count_1] += -1 
                        end  
                        end
                    end    
                end                      ##여기서 tuple의 원소를 보고 n_boundary에 update한다. i,j,k가 뭔지는 중요하지않다.세로열만 중요.
                global count_1 += 1      ##count_1에 +1 해준다.##
            end         
        end
    end 
    global count_1 = 1
    return n_boundary                      
end
###############################################################
#############2dim Bounmdary map################################
count_1 = 1
function bound_matrix_2(m,s,t)
    n_boundary = zeros(Int64, m^1, m^2)
    for i = 1:m
        for j = 1:m
            abc_ = []
            abc_ = (j)                           ## 1 ## 
            abc_ = [abc_; (-R_1(i,j,m,s,t))]     ## 2 ##   
            abc_ = [abc_; (-R_2(i,j,m,s,t))]     ## 3 ##   
            abc_ = [abc_; (i)]                   ## 4 ## 
                

            for i_1 = 1:m
                for k_1 = 1:4
                if (i_1) == abc_[k_1]
                    n_boundary[i_1, count_1] += 1
                elseif ((-1)*(i_1)) == abc_[k_1]   
                    n_boundary[i_1, count_1] += -1 
                end
                end    
            end                      ##여기서 tuple의 원소를 보고 n_boundary에 update한다. i,j,k가 뭔지는 중요하지않다.세로열만 중요.

            global count_1 += 1      ##count_1에 +1 해준다.##        
        end
    end 
    global count_1 = 1
    return n_boundary                      
end
################################################################