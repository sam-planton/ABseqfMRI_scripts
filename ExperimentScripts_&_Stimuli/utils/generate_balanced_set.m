function permut_Seq = generate_balanced_set(n)
%GENERATE_BALANCED_SET This function returns a set of indices from 1 to n
% such that the two halves contain 1 version of each sequence
   
    if mod(n,2) ~= 0
        error('n need to be an even number')
    end
    
    function permut_Seq = generate_permutation()
        
        coins = ones(1,n/2);
        list1 = randperm(n/2);
        list2 = randperm(n/2);
        permut_Seq = [list1,list2];
        
        for k = 1:(n/2)
            % toss a coin
            plusX_ornot = round(rand(1));
            coins(permut_Seq(k)) = coins(permut_Seq(k)) - plusX_ornot;
            permut_Seq(k) = plusX_ornot*(n/2) + permut_Seq(k);
        end
        
        for k = (n/2+1):n
            permut_Seq(k) = coins(permut_Seq(k))*(n/2) + permut_Seq(k);
        end
    end

    permut_Seq = generate_permutation();

    while mod(permut_Seq((n/2)),(n/2)) == mod(permut_Seq((n/2+1)),(n/2))
        permut_Seq = generate_permutation();
    end

end