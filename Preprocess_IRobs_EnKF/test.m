figure;
m_proj('ortho','lat',48','long',-123');
m_coast('patch','r');
m_grid('linest','-','xticklabels',[],'yticklabels',[]);

patch(.55*[-1 1 1 -1],.25*[-1 -1 1 1]-.55,'w'); 
text(0,-.55,'M\_Map','fontsize',25,'color','b',...
    'verticalalignment','middle','horizontalalignment','center');
