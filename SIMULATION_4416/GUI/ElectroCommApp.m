function ElectroCommApp
% ElectroCommApp
% ElectroComm Visualizer - a standalone version that DOES NOT require the
% Communications Toolbox. Includes: Logic Simulator, Modulation Visualizer,
% Curve Tracer, Export.

    fig = uifigure('Name','ElectroComm Visualizer','Position',[50 50 1200 700]);
    tg = uitabgroup(fig);

    %% --- TAB 1: Logic Simulator ---
    tab1 = uitab(tg,'Title','Logic Simulator');
    ddGate = uidropdown(tab1,'Items',{'AND','OR','NOT','XOR','NAND','NOR'},...
        'Position',[20 630 120 25]);
    efInput = uieditfield(tab1,'text','Position',[160 630 200 25],'Value','101');
    btnSim = uibutton(tab1,'push','Text','Simulate','Position',[380 630 80 25]);
    btnClear1 = uibutton(tab1,'push','Text','Clear','Position',[470 630 80 25]);
    tbl = uitable(tab1,'Position',[20 500 360 100],'ColumnName',{'Input','Output'});
    axLogic = uiaxes(tab1,'Position',[20 50 550 400]); title(axLogic,'Timing Diagram');

    btnSim.ButtonPushedFcn = @(~,~)logicSim();
    btnClear1.ButtonPushedFcn = @(~,~)(cla(axLogic); tbl.Data = {});

    function logicSim
        s = char(efInput.Value);
        s(s==' ') = [];                 % remove spaces
        if any(~ismember(s, ['0','1']))
            uialert(fig,'Input must be a binary string like "1010".','Bad Input');
            return
        end
        binInputs = double(s - '0');
        N = length(binInputs);
        gate = ddGate.Value;

        % compute output (vector for plotting & string for table)
        switch gate
            case 'NOT'
                outVec = double(~binInputs);    % bitwise NOT
            case 'AND'
                v = prod(binInputs);
                outVec = repmat(v,1,N);
            case 'OR'
                v = max(binInputs);
                outVec = repmat(v,1,N);
            case 'XOR'
                v = mod(sum(binInputs),2);
                outVec = repmat(v,1,N);
            case 'NAND'
                v = double(~prod(binInputs));
                outVec = repmat(v,1,N);
            case 'NOR'
                v = double(~max(binInputs));
                outVec = repmat(v,1,N);
            otherwise
                outVec = zeros(1,N);
        end

        tbl.Data = {s, char('0'+outVec)}; % show output as 0/1 string

        cla(axLogic);
        t = 0:N; % one extra to use stairs easily
        stairs(axLogic,t,[binInputs binInputs(end)],'LineWidth',2); hold(axLogic,'on');
        stairs(axLogic,t,[outVec outVec(end)],'r','LineWidth',2);
        hold(axLogic,'off');
        xlabel(axLogic,'Bit index'); ylim(axLogic,[-0.5 1.5]);
        legend(axLogic,{'Input','Output'});
        grid(axLogic,'on');
    end

    %% --- TAB 2: Modulation Visualizer ---
    tab2 = uitab(tg,'Title','Modulation Visualizer');
    efBits = uieditfield(tab2,'text','Position',[20 630 150 25],'Value','10110');
    ddMod = uidropdown(tab2,'Items',{'ASK','FSK','BPSK','QPSK','16QAM'},...
        'Position',[190 630 120 25]);
    btnShow = uibutton(tab2,'push','Text','Show Signal','Position',[320 630 100 25]);
    btnNoise = uibutton(tab2,'push','Text','Add Noise + BER','Position',[430 630 140 25]);
    axBin = uiaxes(tab2,'Position',[20 350 350 250]); title(axBin,'Binary Input');
    axSig = uiaxes(tab2,'Position',[400 350 350 250]); title(axSig,'Modulated Signal');
    axConst = uiaxes(tab2,'Position',[20 50 350 250]); title(axConst,'Constellation');
    axEye = uiaxes(tab2,'Position',[400 50 350 250]); title(axEye,'Eye Diagram');
    lblBER = uilabel(tab2,'Text','BER: -','Position',[780 630 150 25]);

    btnShow.ButtonPushedFcn = @(~,~)showSignal();
    btnNoise.ButtonPushedFcn = @(~,~)addNoise();

    % shared state for last modulated signal (so Add Noise can use it)
    lastSym = [];
    lastTx = [];
    lastSps = 8;

    function showSignal
        s = char(efBits.Value); s(s==' ')=[];
        if isempty(s) || any(~ismember(s,['0','1']))
            uialert(fig,'Bits must be a nonempty binary string like "101010".','Bad input'); return
        end
        bits = double(s - '0');
        scheme = ddMod.Value;
        sps = 8;
        [sym, tx] = modulate(bits, scheme, sps);

        lastSym = sym; lastTx = tx; lastSps = sps;

        % Plot binary input (one step per bit)
        cla(axBin);
        stairs(axBin,0:length(bits),[bits bits(end)],'LineWidth',2); ylim(axBin,[-0.5 1.5]); grid(axBin,'on');
        xlabel(axBin,'bit index');

        % Plot modulated signal (real part)
        cla(axSig);
        plot(axSig, (0:length(tx)-1)/sps, real(tx),'LineWidth',1.2);
        xlabel(axSig,'time (symbols)'); grid(axSig,'on');

        % Constellation (use symbol-level points)
        cla(axConst);
        scatter(axConst, real(sym), imag(sym), 36, 'filled');
        xlabel(axConst,'In-phase'); ylabel(axConst,'Quadrature'); grid(axConst,'on');
        title(axConst,'Constellation');

        % Eye diagram (real part)
        cla(axEye);
        drawEye(axEye, tx, sps);
    end

    function addNoise
        s = char(efBits.Value); s(s==' ')=[];
        bits = double(s - '0');
        scheme = ddMod.Value;
        if isempty(lastTx)
            % if user hasn't generated signal yet, create it now
            [~,tx] = modulate(bits, scheme, lastSps);
        else
            tx = lastTx; sym = lastSym;
        end
        SNRdB = 10; % default SNR
        rx = add_awgn(tx, SNRdB);

        % re-sample at symbol centers to get symbol-level rx
        sps = lastSps;
        sampleIdx = ceil(sps/2):sps:length(rx);
        rxSym = rx(sampleIdx);

        rxBits = demodulate(rxSym, scheme, length(bits));
        [~, ber] = biterr(bits(:), rxBits(:));
        lblBER.Text = sprintf('BER: %.5f (SNR=%d dB)', ber, SNRdB);

        % plot noisy constellation
        cla(axConst);
        scatter(axConst, real(rxSym), imag(rxSym), 28, 'o', 'filled');
        xlabel(axConst,'In-phase'); ylabel(axConst,'Quadrature'); grid(axConst,'on');
        title(axConst,'Constellation with Noise');

        % re-plot eye of noisy signal (real part)
        cla(axEye);
        drawEye(axEye, rx, sps);
    end

    %% --- TAB 3: Curve Tracer ---
    tab3 = uitab(tg,'Title','Curve Tracer');
    ddComp = uidropdown(tab3,'Items',{'Silicon Diode','Zener Diode','NPN BJT','NMOS'},...
        'Position',[20 630 150 25]);
    sldParam = uislider(tab3,'Position',[200 640 200 3],'Limits',[1 100],'Value',10);
    btnPlot = uibutton(tab3,'push','Text','Plot Curve','Position',[420 630 100 25]);
    btnClear2 = uibutton(tab3,'push','Text','Clear','Position',[540 630 80 25]);
    axCurve = uiaxes(tab3,'Position',[20 50 700 550]); title(axCurve,'IV Characteristics');

    btnPlot.ButtonPushedFcn = @(~,~)plotCurve();
    btnClear2.ButtonPushedFcn = @(~,~)cla(axCurve);

    function plotCurve
        comp = ddComp.Value; V = linspace(0,10,200);
        cla(axCurve);
        switch comp
            case 'Silicon Diode'
                Is=1e-12; Vt=0.025; I=Is*(exp(V./Vt)-1);
                plot(axCurve,V,I,'LineWidth',2);
                xlabel(axCurve,'V (V)'); ylabel(axCurve,'I (A)');
            case 'Zener Diode'
                Vz=5; Is=1e-12; Vt=0.025; I=Is*(exp(V./Vt)-1);
                I(V<-Vz)=-(V(V<-Vz)+Vz)/10;
                plot(axCurve,V,I,'LineWidth',2);
                xlabel(axCurve,'V (V)'); ylabel(axCurve,'I (A)');
            case 'NPN BJT'
                Beta=round(sldParam.Value);
                Ib=0:10e-6:100e-6; Vce=linspace(0,10,80);
                [IbGrid,VceGrid]=meshgrid(Ib,Vce);
                Ic=Beta.*IbGrid.*(1-exp(-VceGrid/0.2));
                surf(axCurve,VceGrid,IbGrid,Ic);
                xlabel(axCurve,'Vce (V)'); ylabel(axCurve,'Ib (A)'); zlabel(axCurve,'Ic (A)');
            case 'NMOS'
                k=0.5e-3; Vth=2; Vgs=round(sldParam.Value);
                Id=zeros(size(V));
                for i=1:length(V)
                    if V(i) < Vgs-Vth
                        Id(i)=k*((Vgs-Vth)*V(i)-0.5*V(i)^2);
                    else
                        Id(i)=0.5*k*(Vgs-Vth)^2;
                    end
                end
                plot(axCurve,V,Id,'LineWidth',2);
                xlabel(axCurve,'V (V)'); ylabel(axCurve,'Id (A)');
        end
        grid(axCurve,'on');
    end

    %% --- TAB 4: Export ---
    tab4 = uitab(tg,'Title','Export');
    btnExport = uibutton(tab4,'push','Text','Export Plots','FontSize',14,...
        'Position',[400 300 200 50]);
    btnExport.ButtonPushedFcn = @(~,~)exportAll();

    function exportAll
        try
            exportgraphics(axLogic,'logic.png');
        catch; end
        try
            exportgraphics(axSig,'modulated.png');
        catch; end
        try
            exportgraphics(axConst,'constellation.png');
        catch; end
        try
            exportgraphics(axCurve,'curve.png');
        catch; end
        uialert(fig,'All plots exported (where present).','Export Done');
    end

    %% -------------------------
    % Helper functions (modulation/demod/noise/eye)
    %% -------------------------

    function [sym, tx] = modulate(bits, scheme, sps)
        % return symbol-level vector sym and upsampled tx (rect pulses)
        switch upper(scheme)
            case 'ASK'
                sym = 2*bits-1;                  % levels {-1, +1}
                sym = sym(:).';                  % row
            case 'FSK'
                % create two frequencies as complex exponentials
                f0 = 0.1; f1 = 0.3; fs = sps;
                tbit = (0:(sps-1))/fs;
                tx = [];
                for b = bits
                    f = f0*(b==0) + f1*(b==1);
                    tx = [tx exp(1j*2*pi*f*tbit)]; %#ok<AGROW>
                end
                % for consistency, build sym as one complex per bit (center sample)
                sym = tx(ceil(sps/2):sps:end);
                return
            case {'BPSK','PSK'}
                sym = (2*bits-1) + 0j;
            case 'QPSK'
                % group bits into pairs (pad if necessary)
                nb = length(bits);
                pad = mod(2-nb,2);
                if pad>0
                    bits = [bits zeros(1,pad)];
                end
                pairs = reshape(bits,2,[]);
                map = [1+1j, -1+1j, 1-1j, -1-1j]/sqrt(2); % 00,01,10,11
                sym = zeros(1,size(pairs,2));
                for k=1:size(pairs,2)
                    idx = pairs(1,k)*2 + pairs(2,k); % 0..3
                    sym(k) = map(idx+1);
                end
            case {'16QAM','16QAM'}
                % group into 4 bits
                nb = length(bits);
                pad = mod(4-nb,4);
                if pad>0, bits = [bits zeros(1,pad)]; end
                nsyms = length(bits)/4;
                sym = zeros(1,nsyms);
                levels = [-3 -1 1 3];
                for k=1:nsyms
                    nib = bits((k-1)*4 + (1:4));
                    val = nib(1)*8 + nib(2)*4 + nib(3)*2 + nib(4);
                    i_idx = floor(val/4)+1;
                    q_idx = mod(val,4)+1;
                    sym(k) = (levels(i_idx) + 1j*levels(q_idx))/sqrt(10); % normalize
                end
            otherwise
                sym = (2*bits-1)+0j;
        end
        % upsample with rectangular pulses
        tx = repelem(sym, sps);
    end

    function bits_est = demodulate(rxSym, scheme, nbits_expected)
        % rxSym is symbol-level complex array (one complex value per symbol)
        switch upper(scheme)
            case 'ASK'
                bits_est = real(rxSym) > 0;
            case {'BPSK','PSK'}
                bits_est = real(rxSym) > 0;
            case 'QPSK'
                % map back by quadrant signs
                bits_est = [];
                for k=1:length(rxSym)
                    r = rxSym(k);
                    b1 = real(r) < 0; % 1 if negative
                    b2 = imag(r) < 0;
                    bits_est = [bits_est, double(~b1), double(~b2)]; %#ok<AGROW>
                end
                bits_est = bits_est(1:nbits_expected);
            case 'FSK'
                % simple demod: compare instantaneous phase derivative (freq)
                % sample is complex exponential at f0 or f1; decide by sign of mean freq
                bits_est = zeros(1,length(rxSym));
                for k=1:length(rxSym)
                    s = rxSym(k);
                    bits_est(k) = double(angle(s) > 0); % crude
                end
                bits_est = bits_est(1:nbits_expected);
            case '16QAM'
                levels = [-3 -1 1 3]/sqrt(10);
                bits_est = [];
                for k=1:length(rxSym)
                    r = rxSym(k);
                    % find nearest I and Q levels
                    [~, idxI] = min(abs(real(r)-levels));
                    [~, idxQ] = min(abs(imag(r)-levels));
                    i_idx = idxI-1; q_idx = idxQ-1; % 0..3
                    val = i_idx*4 + q_idx;
                    nib = [bitget(val,4) bitget(val,3) bitget(val,2) bitget(val,1)]; % MSB..LSB
                    bits_est = [bits_est nib]; %#ok<AGROW>
                end
                bits_est = bits_est(1:nbits_expected);
            otherwise
                bits_est = real(rxSym) > 0;
        end
        bits_est = bits_est(:);
    end

    function rx = add_awgn(tx, SNRdB)
        % simple AWGN for both real/complex signals
        P = mean(abs(tx(:)).^2);
        N0 = P / (10^(SNRdB/10));
        sigma = sqrt(N0/2);
        rx = tx + sigma*(randn(size(tx)) + 1j*randn(size(tx)));
    end

    function drawEye(ax, tx, sps)
        % draw eye diagram on ax using real(tx)
        sig = real(tx);
        L = length(sig);
        spanSym = 2; % 2-symbol eye
        samplesPerEye = spanSym * sps;
        nEyes = floor(L / samplesPerEye);
        if nEyes < 1
            plot(ax, sig); return
        end
        trimmed = sig(1:nEyes*samplesPerEye);
        M = reshape(trimmed, samplesPerEye, nEyes);
        t = (0:samplesPerEye-1)/sps;
        for k=1:nEyes
            plot(ax, t, M(:,k)); hold(ax,'on');
        end
        hold(ax,'off');
        xlabel(ax,'symbol time'); grid(ax,'on');
    end

end
