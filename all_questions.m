clc; clear; close all;
E0 = -5.5; 
Ef = -5;
gam1 = 0.1 ; 
gam2 = 0.1;
U = 0 ; %Constants (all MKS, except energy which is in eV) 
hbar = 1.06e-34;q = 1.6e-19;IE = (2*q*q)/hbar;
kT = .025; % Bias (calculate 101 voltage points in [-4 4] range) 
nV = 101;VV = linspace(-4,4,nV);dV = VV(2)-VV(1); N0 = 2/ (1 + exp((E0-Ef)/kT) ) ; 
for iV = 1:nV% Voltage loop 
    UU = 0;dU = 1;

    V = VV(iV);mu1 = Ef-(V/2);mu2 = Ef+(V/2); 

    while dU>1e-6%SCF 
        E = E0+UU;

        f1 = 1/(1+exp((E-mu1)/kT));f2 = 1/(1+exp((E-mu2)/kT)); NN = 2* ( (gam1*f1) + (gam2*f2))/(gam1 + gam2) ; % Charge 
        Uold = UU;UU = Uold+(.05*((U*(NN-N0))-Uold)); dU = abs (UU-Uold) ; [V UU dU] ;

    end

    curr = IE*gam1*gam2 *(f2-f1)/(gam1+gam2) ; 
    II(iV) = curr ;
    N(iV) = NN; 
end

[V NN] ;
% voltU0 = figure;



%%
figs = gobjects(6,1);
figs(1,1) = figure('visible','off');

% saveas(figs(1,1), 'newout', 'fig');
%%

G = diff(II)/dV;GG = [G(1) G] ; % Conductance 
voltU0plot = plot(VV,II*10^6,'k');% Plot I-V
xlim([-2 2])

figs(2,1) = figure('visible','off');
% condU0 = figure;
condU0plot = plot(VV,GG*10^6, 'k');
xlim([-2 2])
% saveas(figs(1,1), 'voltU0', 'fig');
%%


U = 1 ; %Constants (all MKS, except energy which is in eV) 
hbar = 1.06e-34;q = 1.6e-19;IE = (2*q*q)/hbar;
kT = .025; % Bias (calculate 101 voltage points in [-4 4] range) 
nV = 101;VV = linspace(-4,4,nV);dV = VV(2)-VV(1); N0 = 2/ (1 + exp((E0-Ef)/kT) ) ; 
for iV = 1:nV% Voltage loop 
    UU = 0;dU = 1;

    V = VV(iV);mu1 = Ef-(V/2);mu2 = Ef+(V/2); 

    while dU>1e-6%SCF 
        E = E0+UU;

        f1 = 1/(1+exp((E-mu1)/kT));f2 = 1/(1+exp((E-mu2)/kT)); NN = 2* ( (gam1*f1) + (gam2*f2))/(gam1 + gam2) ; % Charge 
        Uold = UU;UU = Uold+(.05*((U*(NN-N0))-Uold)); dU = abs (UU-Uold) ; [V UU dU] ;

    end

    curr = IE*gam1*gam2 *(f2-f1)/(gam1+gam2) ; 
    II(iV) = curr ;
    N(iV) = NN; 
end

[V NN] ;
% voltU1 = figure;
figs(3,1) = figure('visible','off');
G = diff(II)/dV;GG = [G(1) G] ; % Conductance 
voltU1plot = plot(VV,II*10^6,'k');% Plot I-V
xlim([-2 2])

figs(4,1) = figure('visible','off');
% condU1 = figure;
condU1plot = plot(VV,GG*10^6, 'k');
xlim([-2 2])
%%
[VV, II, GG] = q9plzwork(0);

[VV2, II2, GG2] = q9plzwork(1);

figs(5,1) = figure('visible','off');
% voltVScurrent = figure;
plot(VV,II*10^6,'G', VV2,II2*10^6,'r')
ylabel('Current(\mu A)')
xlabel('Voltage (V)')
title('Voltage vs Current')
legend('U = 0', 'U = 1')
xlim([-2 2])

% voltVSconductance = figure;
figs(6,1) = figure('visible','off');
voltPlot = plot(VV,GG*10^6,'G',VV2,GG2*10^6, 'r');
ylabel('dI/dV (\mu A/V)')
xlabel('Voltage (V)')
title('Voltage vs Conductance')
legend('U = 0', 'U = 1')
xlim([-2 2])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 10

E0 = -5.5; 
Ef = -5;
gam1 = 0.1 ; 
gam2 = 0.1;
U = 0 ; %Constants (all MKS, except energy which is in eV) 
hbar = 1.06e-34;q = 1.6e-19;IE = (2*q*q)/hbar;
kT = .025; % Bias (calculate 101 voltage points in [-4 4] range) 
nV = 101;VV = linspace(-4,4,nV);dV = VV(2)-VV(1); N0 = 2/ (1 + exp((E0-Ef)/kT) ) ; 
for iV = 1:nV% Voltage loop 
    UU = 0;dU = 1;

    V = VV(iV);mu1 = Ef-(V/2);mu2 = Ef+(V/2); 

    while dU>1e-6%SCF 
        E = E0+UU;

        f1 = 1/(1+exp((E-mu1)/kT));f2 = 1/(1+exp((E-mu2)/kT)); NN = 2* ( (gam1*f1) + (gam2*f2))/(gam1 + gam2) ; % Charge 
        Uold = UU;UU = Uold+(.05*((U*(NN-N0))-Uold)); dU = abs (UU-Uold) ; [V UU dU] ;

    end

    curr = IE*gam1*gam2 *(f2-f1)/(gam1+gam2) ; 
    II(iV) = curr ;
    N(iV) = NN; 
end

[V NN] ;

figs10 = gobjects(7,1);
%%
figs10(1,1) = figure('visible','off');

G = diff(II)/dV;GG = [G(1) G] ; % Conductance 
voltU0plot = plot(VV,II*10^6,'k');% Plot I-V
xlim([-2 2])
%%
figs10(2,1) = figure('visible','off');
condU0plot = plot(VV,GG*10^6, 'k');
xlim([-2 2])

%%


U = 1 ; %Constants (all MKS, except energy which is in eV) 
hbar = 1.06e-34;q = 1.6e-19;IE = (2*q*q)/hbar;
kT = .025; % Bias (calculate 101 voltage points in [-4 4] range) 
nV = 101;VV = linspace(-4,4,nV);dV = VV(2)-VV(1); N0 = 2/ (1 + exp((E0-Ef)/kT) ) ; 
for iV = 1:nV% Voltage loop 
    UU = 0;dU = 1;

    V = VV(iV);mu1 = Ef-(V/2);mu2 = Ef+(V/2); 

    while dU>1e-6%SCF 
        E = E0+UU;

        f1 = 1/(1+exp((E-mu1)/kT));f2 = 1/(1+exp((E-mu2)/kT)); NN = 2* ( (gam1*f1) + (gam2*f2))/(gam1 + gam2) ; % Charge 
        Uold = UU;UU = Uold+(.05*((U*(NN-N0))-Uold)); dU = abs (UU-Uold) ; [V UU dU] ;

    end

    curr = IE*gam1*gam2 *(f2-f1)/(gam1+gam2) ; 
    II(iV) = curr ;
    N(iV) = NN; 
end

[V NN] ;
%%
figs10(3,1) = figure('visible','off');
G = diff(II)/dV;GG = [G(1) G] ; % Conductance 
voltU1plot = plot(VV,II*10^6,'k');% Plot I-V
xlim([-2 2])
%%
figs10(4,1) = figure('visible','off');
condU1plot = plot(VV,GG*10^6, 'k');
xlim([-2 2])
%%
[VV, II, GG] = q9plzwork(0);

[VV2, II2, GG2] = q9plzwork(1);

%%
figs10(5,1) = figure('visible','off');

plot(VV,II*10^6,'G', VV2,II2*10^6,'r')
ylabel('Current(\mu A)')
xlabel('Voltage (V)')
title('Voltage vs Current')
legend('U = 0', 'U = 1')
xlim([-2 2])
%%
figs10(6,1) = figure('visible','off');
voltPlot = plot(VV,GG*10^6,'G',VV2,GG2*10^6, 'r');
ylabel('dI/dV (\mu A/V)')
xlabel('Voltage (V)')
title('Voltage vs Conductance')
legend('U = 0', 'U = 1')
xlim([-2 2])
Ef1 = -2.5; Ef2 = -3.5; Ef3 = -5; 
[VV3, II3, GG3] = q10luckyNo2(Ef1);
[VV4, II4, GG4] = q10luckyNo2(Ef2);
[VV5, II5, GG5] = q10luckyNo2(Ef3);
%%
figs10(7,1) = figure('visible','off');
voltPlot = plot(VV3,II3*10^6,'G',VV4,II4*10^6, 'r', VV5,II5*10^6,'k');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Document generator
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import mlreportgen.report.*                                                 %Import Report API and DOM API packages
import mlreportgen.dom.*

disp("Generating Project 3 Report...")

doctype = "pdf";   %Specify document type for report

styles = containers.Map;
styles("baseHeadingPara") = {Color("darkblue"),FontFamily("Arial")};
styles("heading2Para") = [styles("baseHeadingPara"),{OutlineLevel(2),...
                          OuterMargin("0in","0in","12pt","5pt"),Italic,Underline,FontSize("14pt")}];

                      
rpt = Report("HW3 Report",doctype);                                %Create Report container
rpt.Layout.PageNumberFormat = "i";

tp = TitlePage;                                                             %Create title page reporter
tp.Title = "HW3 Report";                             %Specify title
tp.Author= "Eddie, Summer, Ani";

add(rpt,tp)                                                                 %Add title page to the report    

heading1Para = Paragraph("Question 9 Graphs");
heading1Para.Style = styles("heading2Para");
append(rpt,heading1Para);

% p = Paragraph("Question 9 Graphs");
% append(rpt,p);
    
fig1 = Figure(figs(1,1));                                %Formating image for report
fig1.Scaling = 'custom';
fig1.SnapshotFormat = 'svg';
img = Image(getSnapshotImage(fig1, rpt));
img.Height = "4in";
img.Width = "4.5in";
append(rpt,img);

fig1 = Figure(figs(2,1));                                %Formating image for report
fig1.Scaling = 'custom';
fig1.SnapshotFormat = 'svg';
img = Image(getSnapshotImage(fig1, rpt));
img.Height = "4in";
img.Width = "4.5in";
append(rpt,img);

fig1 = Figure(figs(3,1));                                %Formating image for report
fig1.Scaling = 'custom';
fig1.SnapshotFormat = 'svg';
img = Image(getSnapshotImage(fig1, rpt));
img.Height = "4in";
img.Width = "4.5in";
append(rpt,img);

fig1 = Figure(figs(4,1));                                %Formating image for report
fig1.Scaling = 'custom';
fig1.SnapshotFormat = 'svg';
img = Image(getSnapshotImage(fig1, rpt));
img.Height = "4in";
img.Width = "4.5in";
append(rpt,img);

fig1 = Figure(figs(5,1));                                %Formating image for report
fig1.Scaling = 'custom';
fig1.SnapshotFormat = 'svg';
img = Image(getSnapshotImage(fig1, rpt));
img.Height = "4in";
img.Width = "4.5in";
append(rpt,img);

fig1 = Figure(figs(6,1));                                %Formating image for report
fig1.Scaling = 'custom';
fig1.SnapshotFormat = 'svg';
img = Image(getSnapshotImage(fig1, rpt));
img.Height = "4in";
img.Width = "4.5in";
append(rpt,img);

p = Paragraph("9.c: The conductance gap originates from the voltage division equation, η, when solving for zero bias. The magnitude is determined by the separation of the LUMO and the HOMO and the Fermi energy.");
append(rpt,p);

p = Paragraph(''); append(rpt,p);


p = Paragraph("9.e: This may have something to do with tau becoming greater than the energy. It could also have something to do with a blockage occurring or something about the excited states.");
append(rpt,p);

p = Paragraph(''); append(rpt,p);

heading1Para = Paragraph("Question 10 Answers");
heading1Para.Style = styles("heading2Para");
append(rpt,heading1Para);

% p = Paragraph("Chapter 10 Questions");
% p.Bold;
% append(rpt,p);

p = Paragraph(''); append(rpt,p);

for ff = 1:length(figs10)
    fig1 = Figure(figs10(ff,1));                                %Formating image for report
    fig1.Scaling = 'custom';
    fig1.SnapshotFormat = 'svg';
    img = Image(getSnapshotImage(fig1, rpt));
    img.Height = "4in";
    img.Width = "4.5in";
    append(rpt,img);
end


p = Paragraph("10.d) If n does not equal 0.5 then Cs does not equal Cd. The contacts are no longer symmetric.");
append(rpt,p);

p = Paragraph(''); append(rpt,p);

p = Paragraph("10.e)The device is now in the region where there is no charging, which means Vds=0 and there is no charge in the LUMO. The chemical potential of the source is less than the LUMO energy. ");
append(rpt,p);




close(rpt)                                                          %Close the report
% docview("Project3_Hackett_Report.docx",...
%     "updatefields","showdocxaspdf","closeapp")                      % View the report and convert to PDF
% rptview(rpt)                                                      % Another way to view the report
disp("DONE.")




%%
function [VV, II, GG] = q9plzwork(U)
    E0 = -5.5; 
    Ef = -5;
    gam1 = 0.1 ; 
    gam2 = 0.1;
%     U = 0 ; %Constants (all MKS, except energy which is in eV) 
    hbar = 1.06e-34;q = 1.6e-19;IE = (2*q*q)/hbar;
    kT = .025; % Bias (calculate 101 voltage points in [-4 4] range) 
    nV = 101;VV = linspace(-4,4,nV);dV = VV(2)-VV(1); N0 = 2/ (1 + exp((E0-Ef)/kT) ) ; 
    for iV = 1:nV% Voltage loop 
        UU = 0;dU = 1;
    
        V = VV(iV);mu1 = Ef-(V/2);mu2 = Ef+(V/2); 
    
        while dU>1e-6%SCF 
            E = E0+UU;
    
            f1 = 1/(1+exp((E-mu1)/kT));f2 = 1/(1+exp((E-mu2)/kT)); NN = 2* ( (gam1*f1) + (gam2*f2))/(gam1 + gam2) ; % Charge 
            Uold = UU;UU = Uold+(.05*((U*(NN-N0))-Uold)); dU = abs (UU-Uold) ; [V UU dU] ;
    
        end
    
        curr = IE*gam1*gam2 *(f2-f1)/(gam1+gam2) ; 
        II(iV) = curr ;
        N(iV) = NN; 
    end
    
    [V NN] ;
%     voltU0 = figure;
    G = diff(II)/dV;GG = [G(1) G] ; % Conductance 
%     voltU0plot = plot(VV,II*10^6,'k');% Plot I-V
end

function [VV, II, GG] = q10luckyNo2(Ef)

%     Ef = -5;
    n = 0.5;
    E0 = -1.5; %[-5.5 -1.5];
    gam1 = 0.1; %[.2 ,2];
    gam2 = 0.1; %[.2 .2];
    
    U = 1*[1 1;1 1];
    
    % Constants (all MKS, except energy which is in eV) 
    hbar = 1.06e-34;
    q = 1.6e-19;
    IE = (2*q*q)/hbar;
    kT = .025; 
    n0 = 2./1 + exp((E0-Ef)./kT);
    
    nV = 101; VV = linspace (-6 , 6 , nV) ; 
    dV = VV(2)-VV(1) ;
    Usc = 0; 
    
    for iV = 1:nV 
        dU = 1;
        UU = 0;
    
        V = VV (iV) ; mu1 = Ef+(V/2);mu2 = Ef-(V/2); 
        while dU>1e-6
            E = E0+UU;
            f1 = 1./(1 + exp( (E-mu1) ./kT) );
            f2 = 1./(1 + exp( (E-mu2)./kT) ) ;
        
            Uold = Usc;
        
            Usc = Uold+(.1*(( (n-n0)*U')-Uold) ) ;
            
%             UU = Uold+(.05*((U*(n-n0))-Uold)); 
            dU = abs(Usc-Uold) ; %[V UU dU] ;
        end
        curr = IE*gam1*gam2 *(f2-f1)/(gam1+gam2) ; 
        II(iV) = curr ;
%         N(iV) = NN; 
    
    
%     h = plot(VV,II) ; % Plot I-V
    end
    G = diff (II)/dV;
    GG = [G(1) G] ; 

end