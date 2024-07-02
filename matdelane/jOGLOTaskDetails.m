function jOGLOTaskDetails(nwb)

    blocks = nwb.intervals.get("omission_glo_passive").vectordata.get("task_block_number").data(:);
    correct = nwb.intervals.get("omission_glo_passive").vectordata.get("correct").data(:);
    omission = nwb.intervals.get("omission_glo_passive").vectordata.get("is_omission").data(:);
    stims = nwb.intervals.get("omission_glo_passive").vectordata.get("stimulus_number").data(:);
    
    b1 = (blocks == 1) & (correct == 1) & (stims == 3);

    b2 = (blocks == 2) & (correct == 1) & (stims == 3);
    b2l = (blocks == 2) & (correct == 1) & (isnan(omission)) & (stims == 3);
    b2o2 = (blocks == 2) & (correct == 1) & (omission == 1) & (stims == 3);
    b2o3 = (blocks == 2) & (correct == 1) & (omission == 1) & (stims == 4);
    b2o4 = (blocks == 2) & (correct == 1) & (omission == 1) & (stims == 5);
    
    b3 = (blocks == 3) & (correct == 1) & (stims == 3);

    b4 = (blocks == 4) & (correct == 1) & (stims == 3);
    b4l = (blocks == 4) & (correct == 1) & (isnan(omission)) & (stims == 3);
    b4o2 = (blocks == 4) & (correct == 1) & (omission == 1) & (stims == 3);
    b4o3 = (blocks == 4) & (correct == 1) & (omission == 1) & (stims == 4);
    b4o4 = (blocks == 4) & (correct == 1) & (omission == 1) & (stims == 5);
    
    b5 = (blocks == 5) & (correct == 1) & (stims == 3);
    b5l = (blocks == 5) & (correct == 1) & (isnan(omission)) & (stims == 3);
    b5o2 = (blocks == 5) & (correct == 1) & (omission == 1) & (stims == 3);
    b5o3 = (blocks == 5) & (correct == 1) & (omission == 1) & (stims == 4);
    b5o4 = (blocks == 5) & (correct == 1) & (omission == 1) & (stims == 5);

    fprintf("> Session %s\n\n->Total correct trials: %d\n", nwb.identifier, sum(b1 + b2 + b3 + b4 + b5));
    fprintf("-> AAAB(1) = %d, XAAAB(2) = %d, BBBA(3) = %d, XBBBA(4) = %d, XRRRR(5) = %d \n", sum(b1), sum(b2), sum(b3), sum(b4), sum(b5));
    fprintf("\n--> (1) AAAB habituation : %d \n", sum(b1));

    fprintf("\n--> (2) AAAB omission. omission chance : %f \n", sum(b2o2 + b2o3 + b2o4)/sum(b2));
    fprintf("\n--> (2) AAAB : %d \n", sum(b2l));
    fprintf("\n--> (2) AXAB : %d \n", sum(b2o2));
    fprintf("\n--> (2) AAXB : %d \n", sum(b2o3));
    fprintf("\n--> (2) AAAX : %d \n", sum(b2o4));

    fprintf("\n--> (3) BBBA habituation : %d \n", sum(b3));

    fprintf("\n--> (4) BBBA omission. omission chance : %f \n", sum(b4o2 + b4o3 + b4o4)/sum(b4));
    fprintf("\n--> (4) BBBA : %d \n", sum(b4l));
    fprintf("\n--> (4) BXBA : %d \n", sum(b4o2));
    fprintf("\n--> (4) BBXA : %d \n", sum(b4o3));
    fprintf("\n--> (4) BBBX : %d \n", sum(b4o4));

    fprintf("\n--> (5) RRRR omission. omission chance : %f \n", sum(b5o2 + b5o3 + b5o4)/sum(b5));
    fprintf("\n--> (5) RRRR : %d \n", sum(b5l));
    fprintf("\n--> (5) RXRR : %d \n", sum(b5o2));
    fprintf("\n--> (5) RRXR : %d \n", sum(b5o3));
    fprintf("\n--> (5) RRRX : %d \n", sum(b5o4));

end