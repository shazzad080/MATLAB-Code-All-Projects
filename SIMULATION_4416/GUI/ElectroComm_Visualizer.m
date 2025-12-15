% ElectroComm_Visualizer.m
% Multi-tab MATLAB UI skeleton integrating:
% - Digital Logic Simulator
% - Digital Modulation Visualizer
% - Interactive Diode & Transistor Curve Tracer
%
% Save this file and run: ElectroComm_Visualizer
function ElectroComm_Visualizer
    % Main figure
    fig = uifigure('Name','ElectroComm Visualizer','Position',[100 100 1100 700]);

    % Tab group
    tg = uitabgroup(fig,'Position',[10 10 1080 680]);
    tab1 = uitab(tg,'Title','Logic Simulator');
    tab2 = uitab(tg,'Title','Modulation Visualizer');
    tab3 = uitab(tg,'Title','Curve Tracer');
    tab4 = uitab(tg,'Title','Export / Help');

    %% Tab 1: Logic Simulator
    % Controls
    uilabel(tab1,'Position',[15 620 120 22],'Text','Gate Type:');
    ddGate = uidropdown(tab1,'Items',{'AND','OR','NOT','XOR','NAND','NOR'},'Position',[90 620 120 22]);

    uilabel(tab1,'Position',[230 620 60 22],'Text','Inputs:');
    txtInputs = uieditfield(tab1,'text','Position',[295 620 150 22],'Value','101');

    btnSimulate = uibutton(tab1,'push','Text','Simulate','Position',[460 620 90 22],...
        'ButtonPushedFcn',@(s,e)simulateLogic(ddGate,txtInputs,axLogic,tableTruth));

    btnClear = uibutton(tab1,'push','Text','Clear','Position',[560 620 60 22],...
        'ButtonPushedFcn',@(s,e)clearLogic(axLogic,tableTruth));

    % Plot axes and truth table
    axLogic = uiaxes(tab1,'Position',[15 300 520 300],'Title','Timing Diagram');
    tableTruth = uitable(tab1,'Position',[550 300 520 300],'Data',cell(1,2),'ColumnName',{'Input','Output'});

    uilabel(tab1,'Position',[15 270 400 22],'Text','Notes: Enter binary string (e.g., 10110). For NOT gate provide single-bit inputs per row.');

    %% Tab 2: Modulation Visualizer
    uilabel(tab2,'Position',[15 620 120 22],'Text','Bitstream:');
    txtBits = uieditfield(tab2,'text','Position',[90 620 240 22],'Value','101101');

    uilabel(tab2,'Position',[350 620 100 22],'Text','Modulation:');
    ddMod = uidropdown(tab2,'Items',{'ASK','FSK','PSK','QPSK','QAM'},'Position',[420 620 120 22]);

    btnShow = uibutton(tab2,'push','Text','Show Signal','Position',[555 620 100 22],...
        'ButtonPushedFcn',@(s,e)showModulation(txtBits,ddMod,axBinary,axMod,axConst,axEye,txtBER));
    btnNoise = uibutton(tab2,'push','Text','Add Noise','Position',[665 620 80 22],...
        'ButtonPushedFcn',@(s,e)addNoiseMod(axMod,txtBER)) ;

    % Axes: binary, modulated, constellation, eye
    axBinary = uiaxes(tab2,'Position',[15 360 520 240],'Title','Binary Waveform');
    axMod = uiaxes(tab2,'Position',[560 360 520 240],'Title','Modulated Signal');
    axConst = uiaxes(tab2,'Position',[15 60 360 280],'Title','Constellation');
    axEye = uiaxes(tab2,'Position',[385 60 320 280],'Title','Eye Diagram');

    uilabel(tab2,'Position',[720 120 60 22],'Text','BER:');
    txtBER = uilabel(tab2,'Position',[790 120 200 22],'Text','N/A');

    %% Tab 3: Curve Tracer
    uilabel(tab3,'Position',[15 640 120 22],'Text','Component:');
    ddComp = uidropdown(tab3,'Items',{'Silicon Diode','Zener Diode','NPN BJT','NMOS'},'Position',[95 640 150 22]);

    uilabel(tab3,'Position',[260 640 90 22],'Text','Parameter 1:');
    sldP1 = uislider(tab3,'Position',[350 650 300 3],'Limits',[0 10],'Value',1);
    lblP1 = uilabel(tab3,'Position',[660 640 120 22],'Text','P1 = 1');
    sldP1.ValueChangedFcn = @(s,e)updateLabel(s,lblP1,'P1');

    uilabel(tab3,'Position',[15 600 90 22],'Text','Parameter 2:');
    sldP2 = uislider(tab3,'Position',[95 610 300 3],'Limits',[0 10],'Value',1);
    lblP2 = uilabel(tab3,'Position',[410 600 120 22],'Text','P2 = 1');
    sldP2.ValueChangedFcn = @(s,e)updateLabel(s,lblP2,'P2');

    btnPlot = uibutton(tab3,'push','Text','Plot Curve','Position',[540 640 100 22],...
        'ButtonPushedFcn',@(s,e)plotCurve(ddComp,sldP1.Value,sldP2.Value,axCurve,txtParams));
    btnClearC = uibutton(tab3,'push','Text','Clear','Position',[650 640 80 22],...
        'ButtonPushedFcn',@(s,e)clearAxes(axCurve));

    axCurve = uiaxes(tab3,'Position',[15 60 720 520],'Title','I-V / Output Characteristics');
    txtParams = uilabel(tab3,'Position',[750 200 300 200],'Text','Parameters:');

    %% Tab 4: Export / Help
    uibutton(tab4,'push','Text','Export All Plots','Position',[20 620 120 22],...
        'ButtonPushedFcn',@(s,e)exportAll(fig));
    helpText = sprintf(['ElectroComm Visualizer\n\nTabs:\n- Logic Simulator: simulate boolean gates and timing\n'... 
        '- Modulation Visualizer: view modulated waveforms and constellation\n'... 
        '- Curve Tracer: sweep and plot component I-V / output curves\n\n'... 
        'Implementation notes:\nPlug your signal processing and device models into the labelled callbacks.']);
    uitextarea(tab4,'Position',[20 20 1040 580],'Value',helpText,'Editable','off');

    %% Helper nested functions (placeholders)
    function simulateLogic(ddGateField,inputField,ax,table)
        binStr = strtrim(inputField.Value);
        if isempty(binStr) || any(~ismember(binStr,['0','1']))
            uialert(fig,'Invalid binary input. Use a string like 10110.','Input error');
            return
        end
        inputs = double(binStr)-'0';
        gate = ddGateField.Value;
        % Simple sample-based simulation. Expand for multi-bit or circuit netlist.
        switch gate
            case 'AND'
                out = all(reshape(inputs,[],1),2); % trivial placeholder
            case 'OR'
                out = any(reshape(inputs,[],1),2);
            case 'NOT'
                out = ~inputs;
            case 'XOR'
                out = xor(inputs(1),inputs(end));
            otherwise
                out = inputs; % placeholder
        end
        % Plot timing (stairs)
        t = 0:1:length(inputs)-1;
        cla(ax);
        stairs(ax,t,inputs,'LineWidth',1.2); hold(ax,'on');
        stairs(ax,t, out,'LineWidth',1.2); hold(ax,'off');
        legend(ax,{'Input','Output'});
        ax.XLim = [0 length(inputs)];
        % Update truth table (simple)
        data = cell(length(inputs),2);
        for k=1:length(inputs)
            data{k,1} = inputs(k);
            data{k,2} = out(min(k,end));
        end
        table.Data = data;
    end

    function clearLogic(ax,table)
        cla(ax);
        table.Data = {};    
    end

    function showModulation(txtBitsField,ddModField,axb,axm,axc,axe,labBER)
        bits = txtBitsField.Value;
        if isempty(bits) || any(~ismember(bits,['0','1']))
            uialert(fig,'Invalid bitstream. Use a string like 101101.','Input error');
            return
        end
        b = double(bits)-'0';
        % Binary waveform
        cla(axb); t = 0:1/100:length(b); bw = repelem(b,100); plot(axb,t(1:length(bw)),bw); axb.Title.String='Binary Waveform';

        % Modulated signal placeholder: show carrier with phase changes for PSK
        cla(axm); fs = 1000; fc = 50; Tsym=1; t2 = 0:1/fs:Tsym*length(b)-1/fs;
        switch ddModField.Value
            case {'PSK','QPSK'}
                % Simple BPSK mapping
                s = 2*b-1; sig = kron(s,ones(1,fs)).*cos(2*pi*fc*t2);
                plot(axm,t2,sig); axm.Title.String='Modulated Signal (placeholder)';
                % Placeholder constellation
                cla(axc); scatter(axc,real([1 -1]),imag([0 0])); axc.Title.String='Constellation (placeholder)';
            case 'ASK'
                s = b; sig = kron(s,ones(1,fs)).*cos(2*pi*fc*t2);
                plot(axm,t2,sig); axm.Title.String='ASK Signal';
                cla(axc); text(0.1,0.1,'No constellation for ASK', 'Parent', axc);
            case 'FSK'
                f1 = fc; f2 = fc*1.5; sig = zeros(size(t2));
                for k=1:length(b)
                    idx = (k-1)*fs + (1:fs);
                    sig(idx) = cos(2*pi*(f1*(~b(k))+f2*(b(k))).*t2(idx));
                end
                plot(axm,t2,sig); axm.Title.String='FSK Signal';
                cla(axc); text(0.1,0.1,'No constellation for FSK', 'Parent', axc);
            case 'QAM'
                % Placeholder QAM points
                sig = kron((2*b-1),ones(1,fs)).*cos(2*pi*fc*t2);
                plot(axm,t2,sig); axm.Title.String='QAM Signal (placeholder)';
                cla(axc); scatter(axc,[-1 1 -1 1],[-1 -1 1 1]); axc.Title.String='Constellation (placeholder)';
            otherwise
                plot(axm,t2,zeros(size(t2))); axm.Title.String='Modulated Signal';
        end
        % Eye diagram placeholder
        cla(axe); plot(axe,reshape(sig(1:floor(end/ Tsym*1)),[],1)); axe.Title.String='Eye Diagram (placeholder)';
        % BER placeholder
        labBER.Text = 'BER: 0.00 (placeholder)';
    end

    function addNoiseMod(axModField,labBER)
        % Placeholder: overlay noise
        ax = axModField;
        lines = findobj(ax,'Type','line');
        if isempty(lines), return; end
        y = lines(1).YData + 0.2*randn(size(lines(1).YData));
        hold(ax,'on');
        plot(ax,lines(1).XData,y); hold(ax,'off');
        labBER.Text = 'BER: 0.05 (placeholder)';
    end

    function updateLabel(s,label,prefix)
        label.Text = sprintf('%s = %.2f',prefix,s.Value);
    end

    function plotCurve(ddCompField,p1,p2,ax,labelOut)
        comp = ddCompField.Value;
        cla(ax);
        switch comp
            case 'Silicon Diode'
                % Shockley equation
                V = linspace(-1,1,400);
                Is = 1e-12*p1; n = 1 + p2/10; Vt = 25.85e-3;
                I = Is*(exp(V./(n*Vt))-1);
                semilogy(ax,V,abs(I)); ax.Title.String='Diode I-V (log scale)';
                labelOut.Text = sprintf('Is=%.2e, n=%.2f',Is,n);
            case 'Zener Diode'
                V = linspace(-10,2,400);
                Iz = p1*1e-3; Vz = -5 - p2;
                I = zeros(size(V));
                I(V<Vz) = Iz*exp((Vz-V(V<Vz))/0.1);
                plot(ax,V,I); ax.Title.String='Zener I-V';
                labelOut.Text = sprintf('Vz=%.2f V, Iz-scale=%.2f mA',Vz,Iz);
            case 'NPN BJT'
                Vce = linspace(0,10,200);
                Ib = 1e-3 * p1; beta = 50 + 100*p2/10;
                Ic = beta*Ib*ones(size(Vce));
                plot(ax,Vce,Ic); ax.Title.String='BJT Output Curves (approx)';
                xlabel(ax,'V_{CE} (V)'); ylabel(ax,'I_C (A)');
                labelOut.Text = sprintf('I_B=%.3fmA, beta=%.1f',Ib*1e3,beta);
            case 'NMOS'
                Vds = linspace(0,10,200);
                Vgs = 2 + p1; Vth = 1 + p2/10; k = 1e-3;
                Id = zeros(size(Vds));
                for i=1:length(Vds)
                    if Vgs <= Vth
                        Id(i)=0;
                    elseif Vds(i) < (Vgs-Vth)
                        Id(i) = k*((Vgs-Vth)*Vds(i) - 0.5*Vds(i)^2);
                    else
                        Id(i) = 0.5*k*(Vgs-Vth)^2;
                    end
                end
                plot(ax,Vds,Id); ax.Title.String='MOSFET I-V (approx)';
                xlabel(ax,'V_{DS} (V)'); ylabel(ax,'I_D (A)');
                labelOut.Text = sprintf('V_G_S=%.2f V, Vth=%.2f V',Vgs,Vth);
            otherwise
                text(ax,0.1,0.5,'Component not implemented');
        end
    end

    function clearAxes(ax)
        cla(ax);
    end

    function exportAll(mainFig)
        % Export each axes to PNG files
        figs = findall(mainFig,'Type','axes');
        for k=1:length(figs)
            try
                fname = sprintf('export_ax_%d.png',k);
                exportgraphics(figs(k),fname,'Resolution',150);
            catch
            end
        end
        uialert(fig,'Export complete (files saved to MATLAB working folder).','Export');
    end
end
