% colormap for plotting infrared BT
function cmap = IR_colormap(c_int)
    cmap(1:5/c_int,1)=0.5:c_int/(5-c_int)*0.5:1;
    cmap(1:5/c_int,2)=0;
    cmap(1:5/c_int,3)=0.5:c_int/(5-c_int)*0.5:1;
    cmap((5:c_int:10)/c_int,1)=1;
    cmap((5:c_int:10)/c_int,2)=0:c_int/5*0.5:0.5;
    cmap((5:c_int:10)/c_int,3)=1;
    cmap((10:c_int:20)/c_int,1)=0.8:-c_int/10*0.8:0;
    cmap((10:c_int:20)/c_int,2)=0.8:-c_int/10*0.8:0;
    cmap((10:c_int:20)/c_int,3)=0.8:-c_int/10*0.8:0;
    cmap((20:c_int:30)/c_int,1)=0:c_int/10:1;
    cmap((20:c_int:30)/c_int,2)=0;
    cmap((20:c_int:30)/c_int,3)=0;
    cmap((30:c_int:40)/c_int,1)=1;
    cmap((30:c_int:40)/c_int,2)=0:c_int/10:1;
    cmap((30:c_int:40)/c_int,3)=0;
    cmap((40:c_int:50)/c_int,1)=1:-c_int/10:0;
    cmap((40:c_int:50)/c_int,2)=1;
    cmap((40:c_int:50)/c_int,3)=0;
    cmap((50:c_int:60)/c_int,1)=0;
    cmap((50:c_int:60)/c_int,2)=1:-c_int/10:0;
    cmap((50:c_int:60)/c_int,3)=0:c_int/10:1;
    cmap((60:c_int:70)/c_int,1)=0;
    cmap((60:c_int:70)/c_int,2)=0:c_int/10:1;
    cmap((60:c_int:70)/c_int,3)=1;
    cmap((70:c_int:140)/c_int,1)=1:-c_int/70:0;
    cmap((70:c_int:140)/c_int,2)=1:-c_int/70:0;
    cmap((70:c_int:140)/c_int,3)=1:-c_int/70:0;
end

