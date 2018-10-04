function [As]=ident_retardation_FDna(W,Ajw,order,MethOpt,PlotOpt)
%Script for the indentification of a parametric model
%of the A(jw) = B(w)/(jw) - A(w) 
%
%Theta_opt = argmin Sum_i |A(jwi)-P(Theta,jwi)/Q(Theta,jwi)|^2
%
% Use: [As]=ident_retardation_FD(W,Ajw,order,MethOpt,PlotOpt)
%
% Output:
% Ks - the indetified transfer function object 
%
% Inputs:
% W,Ajw - are the vectors of frequencuies and the complex 
%        frequency response A(jw) = B(w)/(jw) - A(w).
% order -order of the parametric model to identify.
% MethOpt - [method,iter] method options for LS freq dom regression: 
%           1-quasilinar, 2-iterative quasilinear
%           3 - nonlinear iterative. (see fit_siso_fresp.m)
% PlotOpt - Ploting options for the indetified model freq response
%           PlotOpt= [PlotFlag,LogScale,wsFactor,wminFactor,wmaxFactor]
%           PlotFlag = 1-plot/0-no plot.
%           LogScale = 1-logaritmic/0-linear.
%           wsFactor = sampling freq as a factor of W(2)-W(1);
%           wminFactor = min freq as a factor of W(1);
%           wmaxFactor = max freq as a factor of W(N);
%
%           Example: PlotOpt=[1,1,0.1,0.5,2], plots in logaritmic scale
%           the sampling freq to compute the frequency response is 
%           ws=0.1*(W(2)-W(1)); wmin = 0.5*W(1), and wmax=2*W(N).
%
%Author: Tristan Perez
%Created: 2007-8-6

%Identification settings:
method     = MethOpt(1);
iter       = MethOpt(2);
PlotFlag   = PlotOpt(1);
LogScale   = PlotOpt(2);
wsFactor   = PlotOpt(3);
wminFactor = PlotOpt(4);
wmaxFactor = PlotOpt(5);
 
%Prepare the data for identification
Scale=max(abs(Ajw));
Ajws = Ajw/Scale;
Fresp= Ajws;%./complex(0,W); %remove the zero at s=0 from the data

%Frequency response fitting
ord_den = order;
ord_num = order;
[Num,Den,Fresp_hat,stable]=fit_siso_fresp(W,Fresp,[],ord_num,ord_den,method,iter);

%Rescale
Num=Scale*Num;
As =tf(Num,Den);

if PlotFlag
    ws=wsFactor*(W(2)-W(1));
    wmax=W(length(W));
    Wext =[wminFactor*W(1):ws:wmaxFactor*wmax];
    Kw_hat=freqresp(Ks,Wext);
    Kw_hat=reshape(Kw_hat,length(Wext),1);
    %log-scale
    figure
    if LogScale,
        subplot(211)
        semilogx(W,20*log10(abs(Ajw)),'o',Wext,20*log10(abs(Kw_hat)),'LineWidth',1.5)
        grid on
        legend('Data',['Estimate, order = ',num2str(order)])
        xlabel('W [rad/s]')
        ylabel('|A(jw)| [dB]')
        title('Frequency response Identification')
        %
        subplot(212)
        semilogx(W,(180/pi)*angle(Ajw),'o',Wext,(180/pi)*angle(Kw_hat),'LineWidth',1.5)
        grid on
        legend('Data','Estimate')
        xlabel('W [rad/s]')
        ylabel('Ang(A(jw)) [deg]')
    else
        subplot(211)
        plot(W,abs(Ajw),'o',Wext,abs(Kw_hat),'LineWidth',1.5)
        grid on
        legend('Data',['Estimate, order = ',num2str(order)])
        xlabel('W [rad/s]')
        ylabel('|A(jw)|')
        title('Frequency response Identification')
        %
        subplot(212)
        plot(W,(180/pi)*angle(Ajw),'o',Wext,(180/pi)*angle(Kw_hat),'LineWidth',1.5)
        grid on
        legend('Data','Estimate')
        xlabel('W [rad/s]')
        ylabel('Ang(A(jw)) [deg]')
    end
end