% This function is is used to allow for the quick selection of the desired
% multi-port vector fitting method.
% Outputs resulting state space system and
% rmse of fit
function [system,rmserr] = Vf_driver_driver(bigH,s,poles,opts,method,tol,min_tol)

if method == "MF"
    [SER,rmserr_mf,bigHfit,opts2] = VFdriver_non_sym(bigH,s,poles,opts);   
    A = SER.A;    B = SER.B;    C = SER.C;    D = SER.D;    E = SER.E;
    system = ss(A,B,C,D);
    rmserr = rmserr_mf;

elseif method == "SIMO"
    [SER,rmserr_simo,bigYfit,opts2]=VFdriver_multi_SIMO(bigH,s,poles,opts);
%     A_c1 = SER_c1.A; B_c1=SER_c1.B; C_c1=SER_c1.C; D_c1 = SER_c1.D;
%     A_c2 = SER_c2.A; B_c2=SER_c2.B; C_c2=SER_c2.C; D_c2 = SER_c2.D;
%     
%     A_simo = blkdiag(A_c1,A_c2);    B_simo = blkdiag(B_c1,B_c2);    C_simo = [C_c1,C_c2];    D_simo = [D_c1,D_c2];
%     system = ss(A_simo,B_simo,C_simo,D_simo);
%     rmserr = rmserr_simo;

    A = SER.A;    B = SER.B;    C = SER.C;    D = SER.D;    E = SER.E;
    system = ss(A,B,C,D);
    rmserr = rmserr_simo;

elseif method == "PCCF"
    
    [SER,rmserr_pccf,bigHfit,opts2]=VFdriver_PCCF(bigH,s,poles,opts,tol,min_tol,1); 
    
    A = SER.A;    B = SER.B;    C = SER.C;    D = SER.D;    E = SER.E;
    system = ss(A,B,C,D);
    rmserr = rmserr_pccf;

elseif method == "orig"

    [SER,rmserr_orig,bigHfit,opts2]=VFdriver(bigH,s,poles,opts); 
    
    A = SER.A;    B = SER.B;    C = SER.C;    D = SER.D;    E = SER.E;
    system = ss(A,B,C,D);
    rmserr = rmserr_orig;   

end


end