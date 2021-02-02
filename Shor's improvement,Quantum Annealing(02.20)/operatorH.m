function H = operatorH(sd)
    if sd == 's'
        H = [[0,0],[1,0]];
    elseif sd =='d'
        H = [[0,0,0,0],[1,0,0,0],[0,0,sqrt(2),0],[0,0,0,sqrt(3)]];   
    end
end