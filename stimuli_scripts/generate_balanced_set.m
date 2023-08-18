function permut_Seq = generate_balanced_set()
%GENERATE_BALANCED_SET This function returns a set of indices from 1 to 14
% such that the two half contain 1 version of each sequence (i and i+7 corresponding to the same sequence)

    function permut_Seq = generate_permutation()
        
        coins = {1,1,1,1,1,1,1};
        list1 = randperm(7);
        list2 = randperm(7);
        permut_Seq = [list1,list2];
        
        for k = 1:7
            % toss a coin
            plus7_ornot = round(rand(1));
            coins{permut_Seq(k)}  = coins{permut_Seq(k)} - plus7_ornot;
            permut_Seq(k) = plus7_ornot*7 + permut_Seq(k);
        end
        
        for k = 8:14
            permut_Seq(k) = coins{permut_Seq(k)}*7 + permut_Seq(k);
        end
    end

permut_Seq = generate_permutation();

while mod(permut_Seq(7),7) == mod(permut_Seq(8),7)
    permut_Seq = generate_permutation();

end

end