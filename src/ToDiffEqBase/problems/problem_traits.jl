isinplace(prob::AbstractSDDEProblem{uType,tType,lType,iip,ND}) where {uType,tType,lType,iip,ND} = iip
is_diagonal_noise(prob::AbstractSDDEProblem{uType,tType,lType,iip,Nothing}) where {uType,tType,lType,iip,ND} = true
