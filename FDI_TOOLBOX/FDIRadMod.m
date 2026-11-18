function [Krad,Ainf_hat] = FDIRadMod(W,A,Ainf,B,FDIopt,Dof)
% Frequency-domain identification of radiation function models of marine
% structures. This function identifies a SISO transfer function
% corresponding to the coupling specified.
%
% Use: [Krad,Ainf_hat]=FDIRadMod(W,A,Ainf,B,FDIopt,Dof)
%
% W      - is the vector of frequencies at which A(w) and B(w) are compued.
% A      - is the vector of frequency dependant added mass coeffiricents A(w),
% Ainf   -  is the infinite frequency added mass coefficient,
% B      - is the vector of potential damping coefficients B(w),
% FDIopt - is a structure with the following fields,
%
%    FDIopt.OrdMax   - Maximum order of transfer function to be considered.
%                      Typical value 20.
%    FDIopt.AinfFlag - if set to 1, the algorithm uses Ainf (3D hydrodynamic data), 
%                      if set to 0, the algorithm estimates Ainf (2D hydrodynamic data);
%    FDIopt.Method   - There are 3 parameter estimation methods (1,2,3). See help
%                      of fit_siso_fresp. Recomended method 2 
%                      (best trade off between accuracy and speed)
%    FDIopt.Iterations - Related to parameter estimation methods 2 and 3. See help of
%                      fit_siso_fresp. Typical value 20.
%    FDIopt.PlotFlag - If set to 1 the function displays the results.
%    FDIopt.LogLin   - logarithmic or linear frequency scale for plotting.
%    FDIopt.wsFactor - Sample faster than the Hydrodynamic code for plotting. 
%                      Typical value 0.1.
%    FDIopt.wminFactor - The minimum frequency to be used in the plot is 
%                      FDIopt.wminFactor*Wmin, where Wmin is the minimum         
%                      frequency of the dataset used for identification.
%                      Typical value 0.1.
%    FDIopt.wmaxFactor - the maximum frequency to be used in the plot is
%                      FDIopt.wmaxFactor*Wmax, where Wmax is the maximum
%                      frequency of the dataset used for identification.
%                      Typical value 5.
%
% Dof [i j] - vector with the coupling to be indentified. i=1,2,..,6. j=1,2,...,6. 
%
% Description:
%
% The function first prompts the user to select the range of freuencies to
% for idnetification. This is done by clicking with the mouse on a
% plot of the added mass and damping. The low-frequency data point should be clicked first, 
% then the high-frequency point, and finally press return. The selected data is
% then re-plotted.
%
% After selecting the freuquency range, the function allows the elimination
% of data wildpoints. This is done by clicking the points that are to be
% eliminated. Click on all the points to be eliminated and then press
% enter. This does not affect the original data, so if a point in deleted
% accidentaly, the function should can be recalled.
%
% After preparing the data, the function estimates a best order
% approximation based on a measure of fitting of added-mass and damping.
% This is based on the coefficient of determnation R2 commonly used in
% statistics. 
%
% The used has the final option of changing the order via a keyboard input,
% or simply to accept the automatically estimated order and exit the function.
%
%
% Created by Tristan Perez (trisan.perez@ntnu.no)
% Date 2009/9/1, Trondheim, Norway.
% Revisions:

%% Prepapre data for identification: Select frequency range and eliminate
[A,B,W]=EditAB(A,B,W);

%% Compute the retardation function Freq response K(jw)
Kw = B+complex(0,W).*(A-Ainf*ones(size(A)));

%% Frequency domain identification
MethOpt =[FDIopt.Method;FDIopt.Iterations]; 
PlotOpt =[0;FDIopt.LogLin;FDIopt.wsFactor; ...
            FDIopt.wminFactor;FDIopt.wmaxFactor];

%% Use Ainf as part of the dataset        
if FDIopt.AinfFlag==1,
    %Auto order estimates
     FDIopt.Ord = 2; %Initial order of the approximation
     [Ks]=ident_retardation_FD(W,Kw,FDIopt.Ord,MethOpt,PlotOpt);
     %
     Kw_hatFD=freqresp(Ks,W);
     Kw_hatFD=reshape(Kw_hatFD,length(W),1);
     Brecfd = real(Kw_hatFD);
     Arecfd = imag(Kw_hatFD)./W+Ainf*ones(size(W));
     
     SSEB = (B-Brecfd)'*(B-Brecfd);
     SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
     R2B = 1 - SSEB/SSTB;
     
     SSEA = (A-Arecfd)'*(A-Arecfd);
     SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
     R2A = 1 - SSEA/SSTA;
     
     r2Thres = 0.99;
     while  (R2B<r2Thres) & (R2A<r2Thres),
     FDIopt.Ord = FDIopt.Ord +1;
     if FDIopt.Ord > FDIopt.OrdMax,
             break
     end
     [Ks]=ident_retardation_FD(W,Kw,FDIopt.Ord,MethOpt,PlotOpt);
     %Compute coeff of determination
     Kw_hatFD=freqresp(Ks,W);
     Kw_hatFD=reshape(Kw_hatFD,length(W),1);
     Brecfd = real(Kw_hatFD);
     Arecfd = imag(Kw_hatFD)./W+Ainf*ones(size(W));
     
     SSEB = (B-Brecfd)'*(B-Brecfd);
     SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
     R2B = 1 - SSEB/SSTB;
     
     SSEA = (A-Arecfd)'*(A-Arecfd);
     SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
     R2A = 1 - SSEA/SSTA;
     end
    
    %Manual order selection
    for l=1:1000,
    %Call radiation identification function
    [Ks]=ident_retardation_FD(W,Kw,FDIopt.Ord,MethOpt,PlotOpt); 
    %Compute retardation frequency response
     ws=FDIopt.wsFactor*(W(2)-W(1));
     wmax=W(length(W));
     wmin=W(1);
     Wext2 =[FDIopt.wminFactor*wmin:ws:FDIopt.wmaxFactor*wmax]'; 
     Kw_hatFD=freqresp(Ks,Wext2);
     Kw_hatFD=reshape(Kw_hatFD,length(Wext2),1);

    %plot retardation Freq. responses 
    Xpl=[0.01 10];
    figure(102) 
    subplot(221)
    semilogx(W,20*log10(abs(Kw)),'o-r',Wext2,20*log10(abs(Kw_hatFD)),'--b','LineWidth',2)
    xlim(Xpl)
    grid on
    legend('K(jw)',['K_{hat}(jw) order ',num2str(FDIopt.Ord)])
    ylabel('|K(jw)|') %N=Kg m/s^2 = (Kg/s) m/s  
    xlabel('Frequency [rad/s]')
    title(['Convolution Model DoF ',num2str(Dof(1)),num2str(Dof(2))])
    subplot(223)
    semilogx(W,(180/pi)*angle(Kw),'o-r',Wext2,(180/pi)*angle(Kw_hatFD),'--b','LineWidth',2)
    xlim(Xpl)
    grid on
    ylabel('Phase K(jw) [deg]') %N=Kg m/s^2 = (Kg/s) m/s  
    xlabel('Frequency [rad/s]')

    %% Reconstruction of added mass and damping
    %Compute coeff of determination
     Kw_hatFDa=freqresp(Ks,W);
     Kw_hatFDa=reshape(Kw_hatFDa,length(W),1);
     Brecfd = real(Kw_hatFDa);
     Arecfd = imag(Kw_hatFDa)./W+Ainf*ones(size(W));
     
     SSEB = (B-Brecfd)'*(B-Brecfd);
     SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
     R2B = 1 - SSEB/SSTB;
     
     SSEA = (A-Arecfd)'*(A-Arecfd);
     SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
     R2A = 1 - SSEA/SSTA;
    
    
     %Construct B and A t plot
     Brecfd = real(Kw_hatFD);
     Arecfd = imag(Kw_hatFD)./Wext2+Ainf*ones(size(Wext2));
     
    figure(102) 
    subplot(224)
    semilogx(W,B,'o-r',Wext2,Brecfd,'-b','LineWidth',2)
    xlim(Xpl)
    grid on
    legend('B',['Best FD ident, order ', num2str(FDIopt.Ord)])
    ylabel('Damping') %N=Kg m/s^2 = (Kg/s) m/s  
    xlabel('Frequency [rad/s]')
    title(['Potential Damping DoF ',num2str(Dof(1)),num2str(Dof(2))])
    %
    subplot(222)
    semilogx(W,A,'o-r',Wext2,Arecfd,[Wext2(1) Wext2(length(Wext2))],[Ainf Ainf],'--b','LineWidth',2)
    xlim(Xpl)
    legend('A',['Aest FD indet, order ', num2str(FDIopt.Ord)],'Ainf')
    ylabel('Added Mass')
    xlabel('Frequency [rad/s]')
    title(['Added Mass DoF ',num2str(Dof(1)),num2str(Dof(2))])
    grid on

    %% Ask for a new iteration
    display('---------------------------------------------------------------------')   
    KK=input(['current order ',num2str((FDIopt.Ord)),', (NEW VALUE (1-20) or RETURN to exit): ']);
    if isempty(KK)
        break
    end
    FDIopt.Ord = KK;
    end
end

Ainf_hat = Ainf;

%% Do not use Ainf as part of the dataset
if FDIopt.AinfFlag==0,
    % Compute Complex Coefficients 
    Ac = (B./complex(0,W)+A);

    % Frequency domain identification
    MethOpt =[FDIopt.Method;FDIopt.Iterations]; 
    PlotOpt =[0;FDIopt.LogLin;FDIopt.wsFactor; ...
                FDIopt.wminFactor;FDIopt.wmaxFactor];

        
    %Auto order estimates
    FDIopt.Ord = 2;
    [As]=ident_retardation_FDna(W,Ac,FDIopt.Ord,MethOpt,PlotOpt);
    %
    Ainf_hat=As.num{1,1}(1);
    P = [As.num{1,1}(:)-As.den{1,1}(:)*Ainf_hat ;0];
    Q = As.den{1,1}(:);
    Ks =tf(P(3:FDIopt.Ord+2,1)',Q'); %Estimated fluid memory model TF.
    %
     Kw_hatFD=freqresp(Ks,W);
     Kw_hatFD=reshape(Kw_hatFD,length(W),1);
     Brecfd = real(Kw_hatFD);
     Arecfd = imag(Kw_hatFD)./W+Ainf*ones(size(W));
     
     SSEB = (B-Brecfd)'*(B-Brecfd);
     SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
     R2B = 1 - SSEB/SSTB;
     
     SSEA = (A-Arecfd)'*(A-Arecfd);
     SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
     R2A = 1 - SSEA/SSTA;
     
     r2Thres = 0.99;
     
     while  (R2B<r2Thres) & (R2A<r2Thres),
         FDIopt.Ord = FDIopt.Ord +1;
         if FDIopt.Ord > FDIopt.OrdMax,
             break
         end
         [As]=ident_retardation_FDna(W,Ac,FDIopt.Ord,MethOpt,PlotOpt);
         
         Ainf_hat=As.num{1,1}(1);
         P = [As.num{1,1}(:)-As.den{1,1}(:)*Ainf_hat ;0];
         Q = As.den{1,1}(:);
         Ks =tf(P(3:FDIopt.Ord+2,1)',Q'); 
         %Compute coeff of determination
         Kw_hatFD=freqresp(Ks,W);
         Kw_hatFD=reshape(Kw_hatFD,length(W),1);
         Brecfd = real(Kw_hatFD);
         %Arecfd = imag(Kw_hatFD)./W+Ainf*ones(size(W));
         Arecfd = imag(Kw_hatFD)./W+Ainf_hat*ones(size(W));
         SSEB = (B-Brecfd)'*(B-Brecfd);
         SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
         R2B = 1 - SSEB/SSTB;

         SSEA = (A-Arecfd)'*(A-Arecfd);
         SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
         R2A = 1 - SSEA/SSTA;
     end
    
    for l=1:1000,

        [As]=ident_retardation_FDna(W,Ac,FDIopt.Ord,MethOpt,PlotOpt);
        %% Compute retardation frequency response
         ws=FDIopt.wsFactor*(W(2)-W(1));
         wmax=W(length(W));
         wmin=W(1);
         Wext2 =[FDIopt.wminFactor*wmin:ws:FDIopt.wmaxFactor*wmax]';
         Ac_hatFD=freqresp(As,Wext2);
         Ac_hatFD=reshape(Ac_hatFD,length(Wext2),1);

        %% Estimated added mass and fluid memory model from As  
        %As is the identified parametric model of Ac
        Ainf_hat=As.num{1,1}(1);
        P = [As.num{1,1}(:)-As.den{1,1}(:)*Ainf_hat ;0];
        Q = As.den{1,1}(:);
        Ks =tf(P(3:FDIopt.Ord+2,1)',Q'); %Estimated fluid memory model TF.

        %Compute retardation frequency response
         ws=FDIopt.wsFactor*(W(2)-W(1));
             wmax=W(length(W));
             wmin=W(1);
             Wext2 =[FDIopt.wminFactor*wmin:ws:FDIopt.wmaxFactor*wmax]';
         Kw_hatFD=freqresp(Ks,Wext2);
         Kw_hatFD=reshape(Kw_hatFD,length(Wext2),1);

        %plot fluid memory model frequency response 
        Xpl=[min(Wext2) max(Wext2)];

        figure(102) 
        subplot(211)
        semilogx(Wext2,20*log10(abs(Kw_hatFD)),'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        legend(['K hat(jw) order ',num2str(FDIopt.Ord)])
        ylabel('|K(jw)|') %N=Kg m/s^2 = (Kg/s) m/s  
        xlabel('Frequency [rad/s]')
        title('Fluid-Memory Model Frequency Response')
        subplot(212)
        semilogx(Wext2,(180/pi)*angle(Kw_hatFD),'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        ylabel('Phase K(jw) [deg]') %N=Kg m/s^2 = (Kg/s) m/s  
        xlabel('Frequency [rad/s]')
        %end

        %% Reconstruction of Added Mass and Damping

        Bhat=real(Kw_hatFD);
        Ahat =imag(Kw_hatFD)./Wext2+Ainf_hat*ones(size(Wext2));

        figure(103) 
        subplot(221)
        semilogx(W,(abs(Ac)),'o-r',Wext2,(abs(Ac_hatFD)),'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        title(['DoF ',num2str(Dof(1)),num2str(Dof(2))])
        legend('Ac(jw)',['Achat(jw) order ',num2str(FDIopt.Ord)])
        ylabel('|Ac(jw)|') %N=Kg m/s^2 = (Kg/s) m/s  
        xlabel('Freq. [rad/s]')
        subplot(223)
        semilogx(W,(180/pi)*angle(Ac),'o-r',Wext2,(180/pi)*angle(Ac_hatFD),'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        ylabel('Phase Ac(jw) [deg]') %N=Kg m/s^2 = (Kg/s) m/s  
        xlabel('Frequency [rad/s]')
        %
        subplot(222)
        semilogx(W,A,'o-r',Wext2,Ahat,'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        legend('A',['Aest'])
        %end
        title(['DoF ',num2str(Dof(1)),num2str(Dof(2))])
        ylabel('Added Mass')
        xlabel('Frequency [rad/s]')
        subplot(224)
        semilogx(W,B,'o-r',Wext2,Bhat,'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        legend('B',['Best'])
        ylabel('Damping')
        xlabel('Frequency [rad/s]')
        %% Ask for a new iteration
        display('---------------------------------------------------------------------')
        if isstable(Ks)==0
            disp('FAILURE UNSTABLE TRANSFER FUNCTION: The identified transfer must be stable')
            disp('Try to reduce the order (e.g. 3,4,5) or abort and use w_max = 2-3 rads/s');
        end
        KK=input(['Current order ',num2str((FDIopt.Ord)),', (NEW VALUE (1-20) or RETURN to exit): ']);
        if isempty(KK)
            break
        end
        FDIopt.Ord=KK;
    end
    
end

%% Return Result
Krad=Ks;


