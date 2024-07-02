function jGLOTaskDetails(nwb)

    a = (nwb.intervals.get("passive_glo").vectordata.get("task_block_number").data(:) == 1) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:) & (nwb.intervals.get("passive_glo").vectordata.get("stimulus_number").data(:) == 5);
    b = (nwb.intervals.get("passive_glo").vectordata.get("task_block_number").data(:) == 2) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:) & (nwb.intervals.get("passive_glo").vectordata.get("stimulus_number").data(:) == 5);
    c = (nwb.intervals.get("passive_glo").vectordata.get("task_block_number").data(:) == 3) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:) & (nwb.intervals.get("passive_glo").vectordata.get("stimulus_number").data(:) == 5);    
    d = (nwb.intervals.get("passive_glo").vectordata.get("task_block_number").data(:) == 4) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:) & (nwb.intervals.get("passive_glo").vectordata.get("stimulus_number").data(:) == 5);
    
    fprintf("> Session %s\n\n->Total correct trials: %d\n", nwb.identifier, sum(a + b + c + d));
    fprintf("-> Habituation(1) = %d, Main(2) = %d, Random control(3) = %d, Sequence control(4) = %d \n", sum(a), sum(b), sum(c), sum(d));
    fprintf("\n--> (1) LO habituation : %d \n", sum(a));

    a = nwb.intervals.get("passive_glo").vectordata.get("go_gloexp").data(:) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:) & (nwb.intervals.get("passive_glo").vectordata.get("task_block_number").data(:) == 2);
    fprintf("\n--> (2) GO main : %d ", sum(a));
    b = nwb.intervals.get("passive_glo").vectordata.get("lo_gloexp").data(:) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:) & (nwb.intervals.get("passive_glo").vectordata.get("task_block_number").data(:) == 2);
    fprintf("\n--> (2) LO main : %d \n", sum(b));

    a = nwb.intervals.get("passive_glo").vectordata.get("go_rndctl").data(:) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:);
    fprintf("\n--> (3) GO random control : %d ", sum(a));
    b = nwb.intervals.get("passive_glo").vectordata.get("lo_rndctl").data(:) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:);
    fprintf("\n--> (3) LO random control : %d ", sum(b));
    
    c = nwb.intervals.get("passive_glo").vectordata.get("igo_rndctl").data(:) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:);
    fprintf("\n--> (3) iGO random control : %d ", sum(c));
    d = nwb.intervals.get("passive_glo").vectordata.get("ilo_rndctl").data(:) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:);
    fprintf("\n--> (3) iLO random control : %d ", sum(d));

    e = nwb.intervals.get("passive_glo").vectordata.get("rndctl").data(:) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:) & (nwb.intervals.get("passive_glo").vectordata.get("stimulus_number").data(:) == 5);
    fprintf("\n--> (3) Other random control : %d \n", sum(e) - sum(a + b + c + d));

    a = nwb.intervals.get("passive_glo").vectordata.get("go_seqctl").data(:) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:);
    fprintf("\n--> (4) GO sequence control : %d ", sum(a));
    b = nwb.intervals.get("passive_glo").vectordata.get("igo_seqctl").data(:) & nwb.intervals.get("passive_glo").vectordata.get("correct").data(:);
    fprintf("\n--> (4) iGO sequence control : %d \n", sum(b));

end