function G = operatorG(sd)
    if sd =='s'
        G = [[0,1],[0,0]];
    elseif sd =='d'
        G = [[0,1,0,0],[0,0,sqrt(2),0],[0,0,0,sqrt(3)],[0,0,0,0]];
    end
end