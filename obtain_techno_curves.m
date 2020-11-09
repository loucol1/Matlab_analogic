function obtain_techno_curves(technos, technos_name)
    gmid = 3:33;
    col = {'b', 'k', 'r', 'c'};
    
    figure;
    
%     gmid = 19;
%     w = 5e-6;
%     Ni = 5;
%     dwi = 0;
%     dli = 0;
%     
%     cg = Ni*(interpn(technos{3}.VDS, technos{3}.VBS, technos{3}.l, technos{3}.GMID, technos{3}.cgg, vdsi, vbsi, Li, gmid)*(w-dwi)*(Li-dli)...
%             + interpn(technos{3}.VDS, technos{3}.VBS, technos{3}.l, technos{3}.GMID, technos{3}.cgso, vdsi, vbsi, Li, gmid)*(w-dwi)...
%             + interpn(technos{3}.VDS, technos{3}.VBS, technos{3}.l, technos{3}.GMID, technos{3}.cgdo, vdsi, vbsi, Li, gmid)*(w-dwi))
         
         
    for iter = 1:length(technos)
		iter
        techno = technos{iter};
        for iter2 = 1:length(gmid)
            in(iter2) = abs(interp1(techno.GMID, techno.IN, gmid(iter2)));
            vdsat(iter2) = abs(interp1(techno.GMID, techno.VDSAT, gmid(iter2)));
            vgs(iter2) = abs(interp1(techno.GMID, techno.VGS, gmid(iter2)));
        end
        
        subplot(2,2,1)
        semilogx(in, gmid, strcat('-', col{iter})); hold on
        
        subplot(2,2,2)
        semilogx(in, vdsat, strcat('-', col{iter})); hold on
        
        subplot(2,2,3)
        semilogx(in, vgs, strcat('-', col{iter})); hold on
        
        subplot(2,2,4)
        plot(gmid, vgs, strcat('-', col{iter})); hold on
        
    end
    
    subplot(2,2,1)
    xlabel('Id/(W/L)');
    ylabel('gm/Id');
    legend(technos_name);
    
    subplot(2,2,2)
    xlabel('Id/(W/L)');
    ylabel('vdsat');
    legend(technos_name);
    
    subplot(2,2,3)
    xlabel('Id/(W/L)');
    ylabel('vgs');
    legend(technos_name);
    
    subplot(2,2,4)
    xlabel('gm/Id');
    ylabel('vgs');
    legend(technos_name);
end


