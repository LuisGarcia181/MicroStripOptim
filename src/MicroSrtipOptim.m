function MicroSrtipOptim
    global app
    x_fig = 50;
    y_fig = 200; 
    ancho_fig = 400;   
    alto_fig = 350;
    if isunix
        alto_fig = 370;
    end

    app.fig = figure('Name','MicroSrtipOptim','Color',[0.9400 0.9400 0.9400],...
                    'MenuBar','none','NumberTitle','off',...
                   'Position',[x_fig y_fig ancho_fig alto_fig]);

    ap = alto_fig-40;%-40
    
    % epsilon_r
    LaTeX_var(app.fig, [10, ap-2, 10, 20],'E')
    LaTeX_sub_var(app.fig, [18, ap-2, 7, 10],'r')  
    app.e_eps_r = uicontrol('Style','edit', 'FontSize',10,...
                'Position',[30, ap, 40, 20],'String','4.4',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    
    % h
    LaTeX_var(app.fig, [10, ap-32, 10, 20],'h')
    app.e_h = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[30, ap-30, 40, 20],'String','1.27',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    LaTeX_label(app.fig, [75, ap-32, 30, 20], 'mm')
    
    % fc     
    LaTeX_var(app.fig, [10, ap-62, 10, 20],'f')
    LaTeX_sub_var(app.fig, [16, ap-62, 7, 10],'c')     
    app.e_fc = uicontrol('Style','edit', 'FontSize',10,...
               'Position',[30, ap-60, 40, 20],'String','2',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    LaTeX_label(app.fig, [75, ap-62, 40, 20], 'GHz')
    
    % Orden del filtro, n 
    LaTeX_var(app.fig, [10, ap-92, 10, 20],'n')
    app.e_n = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[30, ap-90, 40, 20],'String','7',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    
    % Z_0
    LaTeX_var(app.fig, [120, ap-02, 12, 20],'Z')
    LaTeX_sub(app.fig, [130, ap-02, 7, 10],'0')
    app.e_Zo = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[145, ap-0, 40, 20],'String','50',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    LaTeX_label(app.fig, [190, ap-02, 50, 20], 'Ohms'); % 

    
    % WL
    LaTeX_var(app.fig, [120, ap-32, 17, 20],'W')
    LaTeX_sub_var(app.fig, [131, ap-32, 9, 10],'L')
    app.e_WL = uicontrol('Style','edit', 'FontSize',10,...
                'Position',[145, ap-30, 40, 20],'String','0.2',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    LaTeX_label(app.fig, [190, ap-32, 30, 20], 'mm')
    
    % WC
    LaTeX_var(app.fig, [120, ap-62, 17, 20],'W')
    LaTeX_sub_var(app.fig, [131, ap-62, 9, 10],'C')
    app.e_WC = uicontrol('Style','edit', 'FontSize',10,...
                 'Position',[145, ap-60, 40, 20],'String','20',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    LaTeX_label(app.fig, [190, ap-62, 30, 20], 'mm')
    
    % Tipe of filter
    LaTeX_label(app.fig, [110, ap-95, 110, 20], 'Select Filter')
    filtros = {'Butterworth - 0.5 dB', 'Butterworth - 3 dB',...
               'Butterworth - MAXIMALLY FLAT', 'Chebyshev - Lar=0.01 dB',...
               'Chebyshev - Lar=0.04321 dB','Chebyshev - Lar=0.1 dB'};
    app.pm_TF = uicontrol('Style','popupmenu', 'FontSize',10,...
                 'Position',[30, ap-120, 200, 22],'String',filtros,...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    c=5;
    
    % WLmin
    LaTeX_var(app.fig, [10+c, ap-160, 17, 20],'W')
    LaTeX_sub_var(app.fig, [21+c, ap-160, 40, 10],'Lmin')
    app.e_WLmin = uicontrol('Style','edit', 'FontSize',10,...
                  'Position',[10+c, ap-183, 55, 20],'String','0.4',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');    
    % WLstep
    LaTeX_var(app.fig, [80+c, ap-160, 17, 20],'W')
    LaTeX_sub_var(app.fig, [91+c, ap-160, 40, 12],'Lstep')    
    app.e_WLstep = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[80+c, ap-183, 55, 20],'String','0.05',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    % WLmax
    LaTeX_var(app.fig, [150+c, ap-160, 17, 20],'W')
    LaTeX_sub_var(app.fig, [161+c, ap-160, 40, 10],'Lmax')   
    app.e_WLmax = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[150+c, ap-183, 55, 20],'String','1',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    % WCmin
    LaTeX_var(app.fig, [10+c, ap-215, 17, 20],'W')
    LaTeX_sub_var(app.fig, [21+c, ap-215, 40, 10],'Cmin')
    app.e_WCmin = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[10+c, ap-238, 55, 20],'String','2',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    % WCstep
    LaTeX_var(app.fig, [80+c, ap-215, 17, 20],'W')
    LaTeX_sub_var(app.fig, [91+c, ap-215, 40, 12],'Cstep')   
    app.e_WCstep = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[80+c, ap-238, 55, 20],'String','1',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif'); % 0.5
    % WCmax
    LaTeX_var(app.fig, [150+c, ap-215, 17, 20],'W')
    LaTeX_sub_var(app.fig, [161+c, ap-215, 40, 10],'Cmax')    
    app.e_WCmax = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[150+c, ap-238, 55, 20],'String','50',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    
    % Min length
    LaTeX_label(app.fig, [10+c, ap-270, 80, 20],'min(')
    LaTeX_var(app.fig, [39+c, ap-270, 40, 20],'TL')
    LaTeX_label(app.fig, [57+c, ap-270, 6, 20],')')
    app.e_minLen = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[10+c, ap-290, 55, 20],'String','',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');    
    % const length
    LaTeX_var(app.fig, [80+c, ap-270, 40, 20],'TL')
    LaTeX_label(app.fig, [102+c, ap-270, 45, 20],'const')
    app.e_cLen = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[80+c, ap-290, 55, 20],'String','60',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif'); % 0.5
    % Max length
    LaTeX_label(app.fig, [150+c, ap-270, 80, 20],'max(')
    LaTeX_var(app.fig, [181+c, ap-270, 40, 20],'TL')
    LaTeX_label(app.fig, [198+c, ap-270, 6, 20],')')    
    app.e_maxLen = uicontrol('Style','edit', 'FontSize',10,...
              'Position',[150+c, ap-290, 55, 20],'String','',...
              'BackgroundColor',[1 1 1],'FontName','Microsoft Sans Serif');
    

    % ---------------------- Plot Scattering --------------------------
    btn_grc = uicontrol('Style','pushbutton','String', 'Scattering',...
                   'Callback',{@boton_graf_scattering},...
                   'FontSize',10,'Position',[250,ap-0, 100, 26],...
                   'FontName','Microsoft Sans Serif');

    % ---------------------- Print Data --------------------------
    btn_pd = uicontrol('Style','pushbutton','String', 'Print Data',...
                   'Callback',{@boton_print_data},...
                   'FontSize',10,'Position',[250,ap-30, 100, 26],...
                   'FontName','Microsoft Sans Serif');
                   
    % ---------------------- Encontrar límites ----------------------------
    btn_flim = uicontrol('Style','pushbutton','String', 'Find limits',...
                   'Callback',{@boton_find_limits},...
                   'FontSize',10,'Position',[250,ap-60, 100, 26],...
                   'FontName','Microsoft Sans Serif');    
               
    % ---------------------- Graficar limites vs Zo -------------------------
    btn_lvzo = uicontrol('Style','pushbutton','String', 'limits vs Zo',...
                   'Callback',{@boton_lim_Zo},...
                   'FontSize',10,'Position',[250,ap-90, 100, 26],...
                   'FontName','Microsoft Sans Serif');
               
    % ---------------------- Graficar largo vs fc -------------------------
    btn_lvfc = uicontrol('Style','pushbutton','String', 'TL vs fc',...
                   'Callback',{@boton_length_fc},...
                   'FontSize',10,'Position',[250,ap-120, 100, 26],...
                   'FontName','Microsoft Sans Serif');
               
    % ---------------------- Simular --------------------------------
    btn_sim = uicontrol('Style','pushbutton','String', 'Simulate',...
                   'Callback', {@boton_simular},...
                   'FontSize',10,'Position',[250,ap-160, 100, 26],...
                   'FontName','Microsoft Sans Serif');
               
    app.t_porcentaje = uicontrol('Style','text','FontSize',10,...
                       'Position',[350, ap-160, 50, 20],'String','',...
                   'FontName','Microsoft Sans Serif');
               
    % ---------------------- Graficar TL -----------------------------
    app.btn_bgtl = uicontrol('Style','pushbutton','String', "TL, L_A, L_A'",...
                   'Callback', {@boton_graf_length}, ...
                   'FontSize',10,'Enable','off','Position',[250,ap-190, 100, 26],...
                   'FontName','Microsoft Sans Serif');
               
    % ---------------------- Graficar error -----------------------------
    app.btn_lconst = uicontrol('Style','pushbutton','String', 'TL const',...
                   'Callback',{@boton_lconst},'Enable','off',...
                   'FontSize',10,'Position',[250,ap-220, 100, 26],...
                   'FontName','Microsoft Sans Serif');
end

function LaTeX_label(parent, position,string)
    if (ispc)
        uicontrol('Parent',parent,'Style','text','FontSize',10,'FontName','Lucida Bright',...
              'HorizontalAlignment', 'left',...
              'Position',position,'String',string);
    else
        uicontrol('Parent',parent,'Style','text','FontSize',10,'FontName','cmr10',...
              'HorizontalAlignment', 'left',...
              'Position',position,'String',string);
    end          
end

function LaTeX_var(parent, position,string)
    if (ispc)
        uicontrol('Parent',parent,'Style','text','FontSize',10,'FontName','Lucida Bright',...
              'HorizontalAlignment', 'left', 'FontAngle', 'italic',...
              'Position',position,'String',string);
    else
        uicontrol('Parent',parent,'Style','text','FontSize',10,'FontName','cmmi10',...
              'HorizontalAlignment', 'left',...
              'Position',position,'String',string);
    end
end

function LaTeX_sub(parent, position,string)
    if (ispc)
        uicontrol('Parent',parent,'Style','text','FontSize',8,'FontName','Lucida Bright',...
              'HorizontalAlignment', 'left',...
              'Position',position,'String',string);
    else
        uicontrol('Parent',parent,'Style','text','FontSize',8,'FontName','cmr10',...
              'HorizontalAlignment', 'left',...
              'Position',position,'String',string);
    end
end

function LaTeX_sub_var(parent, position,string)
    if (ispc)
        uicontrol('Parent',parent,'Style','text','FontSize',8,'FontName','Lucida Bright',...
              'HorizontalAlignment', 'left', 'FontAngle', 'italic',...
              'Position',position,'String',string);
    else
        uicontrol('Parent',parent,'Style','text','FontSize',8,'FontName','cmmi10',...
              'HorizontalAlignment', 'left',...
              'Position',position,'String',string);
    end
end

function [Zc, l_g, e_re] = Z_c(W, fc, h, e_r) % Calcula impedancia y longitud de onda
    eta = 376.73;% 120*pi; % impedancia característica del vacío
    
    u = W/h;
    a = 1 + (1/49)*log((u^4+(u/52)^2)/(u^4+0.432)) + (1/18.7)*log(1 + (u/18.1)^3);
    b = 0.564*((e_r-0.9)/(e_r + 3))^0.053;
    
    e_re = ((e_r+1)/2) + ((e_r-1)/2)*(1 + (10/u))^(-a*b);
    F = 6 + (2*pi - 6)*exp(-(30.666/u)^0.7528);
    Zc = (eta/(2*pi*sqrt(e_re)))*log(F/u + sqrt(1 + (2/u)^2));
    l_g = 300/(fc*sqrt(e_re));
end

function [g] = param_g(N)
    %------------------------- Butterworth - 0.5 dB -------------------------- %
    % g1   g2     g3     g4     g5     g6     g7     g8     g9     g10
    G(:,:,1) = [
    0.6986 1.0000 0      0      0      0      0      0      0      0; ... %N=1
    1.4029 0.7071 1.9841 0      0      0      0      0      0      0; ... %N=2
    1.5963 1.0967 1.5963 1.0000 0      0      0      0      0      0; ... %N=3
    1.6703 1.1926 2.3661 0.8419 1.9841 0      0      0      0      0; ... %N=4
    1.7058 1.2296 2.5408 1.2296 1.7058 1.0000 0      0      0      0; ... %N=5
    1.7254 1.2479 2.6064 1.3137 2.4758 0.8696 1.9841 0      0      0; ... %N=6
    1.7372 1.2583 2.6381 1.3444 2.6381 1.2483 1.7372 1.0000 0      0; ... %N=7
    1.7451 1.2647 2.6564 1.3590 2.6964 1.3389 2.5093 0.8796 1.9841 0; ... %N=8
    1.7504 1.2690 2.6678 1.3673 2.7239 1.3673 2.6678 1.2690 1.7504 1.0000]; %N=9

    %-------------------------- Butterworth - 3 dB --------------------------- %
    % g1   g2     g3     g4     g5     g6     g7     g8     g9     g10
    G(:,:,2) = [
    1.9953 1.0000 0      0      0      0      0      0      0      0; ... %N=1
    3.1013 0.5339 5.8095 0      0      0      0      0      0      0; ... %N=2
    3.3487 0.7117 3.3487 1.0000 0      0      0      0      0      0; ... %N=3
    3.4389 0.7483 4.3471 0.5920 5.8095 0      0      0      0      0; ... %N=4
    3.4817 0.7618 4.5381 0.7618 3.4817 1.0000 0      0      0      0; ... %N=5
    3.5045 0.7685 4.6061 0.7929 4.4641 0.6933 5.8095 0      0      0; ... %N=6
    3.5182 0.7723 4.6386 0.8039 4.6386 0.7723 3.5182 1.0000 0      0; ... %N=7
    3.5277 0.7745 4.6575 0.8089 4.6990 0.8018 4.4990 0.6073 5.8095 0; ... %N=8
    3.5340 0.7760 4.6692 0.8118 4.7272 0.8118 4.6692 0.7760 3.5340 1.0000]; %N=9

    %---------------------- Butterworth - MAXIMALLY FLAT --------------------- %
    % g1 g2 g3 g4 g5 g6 g7 g8 g9 g10
    G(:,:,3) = [
    2.0000 1.0000 0 0 0 0 0 0 0 0; ... %N=1
    1.4142 1.4142 1.0000 0 0 0 0 0 0 0; ... %N=2
    1.0000 2.0000 1.0000 1.0000 0 0 0 0 0 0; ... %N=3
    0.7654 1.8478 1.8478 0.7654 1.0000 0 0 0 0 0; ... %N=4
    0.6180 1.6180 2.0000 1.6180 0.6180 1.0000 0 0 0 0; ... %N=5
    0.5176 1.4142 1.9318 1.9318 1.4142 0.5176 1.0000 0 0 0; ... %N=6
    0.4450 1.2470 1.8019 2.0000 1.8019 1.2470 0.4450 1.0000 0 0; ... %N=7
    0.3902 1.1111 1.6629 1.9615 1.9615 1.6629 1.1111 0.3902 1.0000 0; ... %N=8
    0.3473 1.0000 1.5321 1.8794 2.0000 1.8794 1.5321 1.0000 0.3473 1.0000]; %N=9
    
    %------------------------ Chebyshev - Lar=0.01 dB ------------------------ %
    % g1   g2     g3     g4     g5     g6     g7     g8     g9     g10
    G(:,:,4) = [
    0.0960 1.0000 0      0      0      0      0      0      0      0; ... %N=1
    0.4489 0.4078 1.1008 0      0      0      0      0      0      0; ... %N=2
    0.6292 0.9703 0.6292 1.0000 0      0      0      0      0      0; ... %N=3
    0.7129 1.2004 1.3213 0.6476 1.1008 0      0      0      0      0; ... %N=4
    0.7563 1.3049 1.5773 1.3049 0.7563 1.0000 0      0      0      0; ... %N=5
    0.7814 1.3600 1.6897 1.5350 1.4970 0.7098 1.1008 0      0      0; ... %N=6
    0.7970 1.3924 1.7481 1.6331 1.7481 1.3924 0.7970 1.0000 0      0; ... %N=7
    0.8073 1.4131 1.7825 1.6833 1.8529 1.6193 1.5555 0.7334 1.1008 0; ... %N=8
    0.8145 1.4271 1.8044 1.7125 1.9058 1.7125 1.8044 1.4271 0.8145 1.0000]; %N=9

    %---------------------- Chebyshev - Lar=0.04321 dB ----------------------- %
    % g1   g2     g3     g4     g5     g6     g7     g8     g9     g10
    G(:,:,5) = [
    0.2000 1.0000 0      0      0      0      0      0      0      0; ... %N=1
    0.6648 0.5445 1.2210 0      0      0      0      0      0      0; ... %N=2
    0.8516 1.1032 0.8516 1.0000 0      0      0      0      0      0; ... %N=3
    0.9314 1.2920 1.5775 0.7628 1.2210 0      0      0      0      0; ... %N=4
    0.9714 1.3721 1.8014 1.3721 0.9714 1.0000 0      0      0      0; ... %N=5
    0.9940 1.4131 1.8933 1.5506 1.7253 0.8141 1.2210 0      0      0; ... %N=6
    1.0080 1.4368 1.9398 1.6220 1.9398 1.4368 1.0080 1.0000 0      0; ... %N=7
    1.0171 1.4518 1.9667 1.6574 2.0237 1.6107 1.7726 0.8330 1.2210 0; ... %N=8
    1.0235 1.4619 1.9837 1.6778 2.0649 1.6778 1.9837 1.4619 1.0235 1.0000]; %N=9
    
    %------------------------ Chebyshev - Lar=0.1 dB ------------------------- %
    % g1 g2 g3 g4 g5 g6 g7 g8 g9 g10
    G(:,:,6) = [
    0.3052 1.0000 0 0 0 0 0 0 0 0; ... %N=1
    0.8431 0.6220 1.3554 0 0 0 0 0 0 0; ... %N=2
    1.0316 1.1474 1.0316 1.0000 0 0 0 0 0 0; ... %N=3
    1.1088 1.3062 1.7704 0.8181 1.3554 0 0 0 0 0; ... %N=4
    1.1468 1.3712 1.9750 1.3712 1.1468 1.0000 0 0 0 0; ... %N=5
    1.1681 1.4040 2.0562 1.5171 1.9029 0.8618 1.3554 0 0 0; ... %N=6
    1.1812 1.4228 2.0967 1.5734 2.0967 1.4228 1.1812 1.0000 0 0; ... %N=7
    1.1898 1.4346 2.1199 1.6010 2.1700 1.5641 1.9445 0.8778 1.3554 0; ... %N=8
    1.1957 1.4426 2.1346 1.6167 2.2054 1.6167 2.1346 1.4426 1.1957 1.0000]; %N=9

    global app
    filtro = get(app.pm_TF,'Value');
    
    g = G(N,1:N,filtro);
end

function [l_L, l_C, l_Zo] = calc_largo(W0, WL, WC, Fc, N, h, e_r)
    fc = Fc*1e9;
    impar = @(x) mod(x,2);
    [Z_0, l_g0] = Z_c(W0, Fc, h, e_r);
    x = 3:6; 
    Lzo = l_g0./(2.^(x-1));
    i = 4;
    while (Lzo(i)<5)&&(i>1)
        i = i-1;
    end
    
    l_Zo = Lzo(i); %l_g0/16; % l_g0/8;  % longitud de la línea de Zo

    [Z_0L, l_gL] = Z_c(WL, Fc, h, e_r);
    [Z_0C, l_gC] = Z_c(WC, Fc, h, e_r);
    
    g = param_g(N);
    g0 = 1;
    Omega_c = 1;
    
    L = zeros(1,N);
    C = L;
    for i = 1:N
        if (impar(i))
            K = (Z_0/g0)*(Omega_c/(2*pi*fc));
            L(i) = K*g(i);
        else
            C(i) = (g0/Z_0)*(Omega_c/(2*pi*fc))*g(i);
        end
    end
  
    wc = 2*pi*fc;
    b_L = 2*pi/l_gL; % beta_L
    b_C = 2*pi/l_gC; % beta_C
   
    l_L = (1/b_L).*asin(wc.*L./Z_0L); % * cambio
    l_C = (1/b_C).*asin(wc.*C.*Z_0C); % * cambio

    if (impar(N))
        m = (N+1)/2;
    else
        m = N/2;
    end
 
    switch N
        case 1
            fun = @(x) (Z_0L*sin(b_L*x(1)) - wc*L(1));
            x0 = l_L(1);
        case 2
            fun = @(x) [Z_0L*sin(b_L*x(1)) + Z_0C*tan(b_C*x(2)/2) - wc*L(1);
                (1/Z_0C)*sin(b_C*x(2)) + (1/Z_0L)*tan(b_L*x(1)/2) - wc*C(2)];
            x0 = l_L(1:N)+l_C(1:N); 
        case 3
            fun = @(x) [Z_0L*sin(b_L*x(1)) + Z_0C*tan(b_C*x(2)/2) - wc*L(1);
                (1/Z_0C)*sin(b_C*x(2)) + (2/Z_0L)*tan(b_L*x(1)/2) - wc*C(2)];
            x0 = l_L(1:m)+l_C(1:m);
        case 4
            fun = @(x) [Z_0L*sin(b_L*x(1)) + Z_0C*tan(b_C*x(2)/2) - wc*L(1);
                (1/Z_0C)*sin(b_C*x(2)) + (1/Z_0L)*(tan(b_L*x(1)/2) + tan(b_L*x(3)/2)) - wc*C(2);
                Z_0L*sin(b_L*x(3)) + Z_0C*(tan(b_C*x(2)/2) + tan(b_C*x(4)/2)) - wc*L(3);
                (1/Z_0C)*sin(b_C*x(4)) + (1/Z_0L)*tan(b_L*x(3)/2) - wc*C(4)];
            x0 = l_L(1:N)+l_C(1:N); 
        case 5
            fun = @(x) [Z_0L*sin(b_L*x(1)) + Z_0C*tan(b_C*x(2)/2) - wc*L(1);
                (1/Z_0C)*sin(b_C*x(2)) + (1/Z_0L)*tan(b_L*x(1)/2) + (1/Z_0L)*tan(b_L*x(3)/2) - wc*C(2);
                Z_0L*sin(b_L*x(3)) + 2*Z_0C*tan(b_C*x(2)/2) - wc*L(3)];
            x0 = l_L(1:m)+l_C(1:m);
        case 6
            fun = @(x) [Z_0L*sin(b_L*x(1)) + Z_0C*tan(b_C*x(2)/2) - wc*L(1);
                (1/Z_0C)*sin(b_C*x(2)) + (1/Z_0L)*(tan(b_L*x(1)/2) + tan(b_L*x(3)/2)) - wc*C(2);
                Z_0L*sin(b_L*x(3)) + Z_0C*(tan(b_C*x(2)/2) + tan(b_C*x(4)/2)) - wc*L(3);
                (1/Z_0C)*sin(b_C*x(4)) + (1/Z_0L)*(tan(b_L*x(3)/2) + tan(b_L*x(5)/2)) - wc*C(4);
                Z_0L*sin(b_L*x(5)) + Z_0C*(tan(b_C*x(4)/2) + tan(b_C*x(6)/2)) - wc*L(5);
                (1/Z_0C)*sin(b_C*x(6)) + (1/Z_0L)*tan(b_L*x(5)/2) - wc*C(6)];
            x0 = l_L(1:N)+l_C(1:N); 
        case 7            
            fun = @(x) [Z_0L*sin(b_L*x(1)) + Z_0C*tan(b_C*x(2)/2) - wc*L(1);
                (1/Z_0C)*sin(b_C*x(2)) + (1/Z_0L)*tan(b_L*x(1)/2) + (1/Z_0L)*tan(b_L*x(3)/2) - wc*C(2);
                Z_0L*sin(b_L*x(3)) + Z_0C*tan(b_C*x(2)/2) + Z_0C*tan(b_C*x(4)/2) - wc*L(3);
                (1/Z_0C)*sin(b_C*x(4)) + (2/Z_0L)*tan(b_L*x(3)/2) - wc*C(4)];
            x0 = l_L(1:m)+l_C(1:m); 
        case 8
            fun = @(x) [Z_0L*sin(b_L*x(1)) + Z_0C*tan(b_C*x(2)/2) - wc*L(1);             % 1
                (1/Z_0C)*sin(b_C*x(2)) + (1/Z_0L)*(tan(b_L*x(1)/2) + tan(b_L*x(3)/2)) - wc*C(2); % 2
                Z_0L*sin(b_L*x(3)) + Z_0C*(tan(b_C*x(2)/2) + tan(b_C*x(4)/2)) - wc*L(3); % 3
                (1/Z_0C)*sin(b_C*x(4)) + (1/Z_0L)*(tan(b_L*x(3)/2) + tan(b_L*x(5)/2)) - wc*C(4); % 4
                Z_0L*sin(b_L*x(5)) + Z_0C*(tan(b_C*x(4)/2) + tan(b_C*x(6)/2)) - wc*L(5); % 5
                (1/Z_0C)*sin(b_C*x(6)) + (1/Z_0L)*(tan(b_L*x(5)/2) + tan(b_L*x(7)/2)) - wc*C(6); % 6
                Z_0L*sin(b_L*x(7)) + Z_0C*(tan(b_C*x(6)/2) + tan(b_C*x(8)/2)) - wc*L(7); % 7
                (1/Z_0C)*sin(b_C*x(8)) + (1/Z_0L)*tan(b_L*x(7)/2) - wc*C(8)];            % 8
            x0 = l_L(1:N)+l_C(1:N); 
        case 9            
            fun = @(x) [Z_0L*sin(b_L*x(1)) + Z_0C*tan(b_C*x(2)/2) - wc*L(1);
                (1/Z_0C)*sin(b_C*x(2)) + (1/Z_0L)*(tan(b_L*x(1)/2) + tan(b_L*x(3)/2)) - wc*C(2);
                Z_0L*sin(b_L*x(3)) + Z_0C*(tan(b_C*x(2)/2) + tan(b_C*x(4)/2)) - wc*L(3);
                (1/Z_0C)*sin(b_C*x(4)) + (1/Z_0L)*(tan(b_L*x(3)/2) + tan(b_L*x(5)/2)) - wc*C(4);
                Z_0L*sin(b_L*x(5)) + 2*Z_0C*tan(b_C*x(4)/2) - wc*L(5)];
            x0 = l_L(1:m)+l_C(1:m); 
    end
    
    options = optimset('Display','off');
    X = fsolve(fun,x0,options);
    
    if (impar(N))
        for j = 1:N
            if (j > m) 
                X(j) = X(2*m-j);
            end
            if (impar(j))
                l_L(j) = X(j);
            else
                l_C(j) = X(j);
            end
        end
    end    
end

function [l_L, l_C, l_Zo, fc] = calc_largo_fc(Zo, W0, WL, WC, Fc, N, h, e_r)
    % esta función calcula los tamaños (largo) efectivos al considerar que
    % en la fc el filtro debe tener una atenuación de 3dB
    [l_L, l_C, l_Zo] = calc_largo(W0, WL, WC, Fc, N, h, e_r); 
    [Y1_dB, S12_fc] = scattering_fc(Zo, WL, WC, N, Fc, h, e_r, l_L, l_C, Fc);
    error_fc = -3 - S12_fc;
    fc = Fc;
    z = 1;
    while (abs(error_fc) > 0.01)&&(z<10)
        fc = fc + 0.02*error_fc;
        [l_L, l_C] = calc_largo(W0, WL, WC, fc, N, h, e_r);
        [Y1_dB, S12_fc] = scattering_fc(Zo, WL, WC, N, fc, h, e_r, l_L, l_C, Fc);
        error_fc = -3 - S12_fc;
        z = z+1;
    end    
end

function [L1, L2, C] = Step_Width(W_1, W_2, Z1, Z2, e_ef1, e_ef2, H)
    if (W_1 > W_2)
        W1 = W_1*1e-3;
        W2 = W_2*1e-3;
        Zc1 = Z1;
        Zc2 = Z2;
        e_eff1 = e_ef1;
        e_eff2 = e_ef2;
        cambio = 0; 
    else
        W1 = W_2*1e-3;
        W2 = W_1*1e-3;
        Zc1 = Z2;
        Zc2 = Z1;
        e_eff1 = e_ef2;
        e_eff2 = e_ef1; 
        cambio = 1;
    end
    h = H*1e3;
    c = 3e8; % the light velocity in free space
    
    Lw1 = Zc1*sqrt(e_eff1)/c; % inductance per unit length of the microstrip 1
    Lw2 = Zc2*sqrt(e_eff2)/c;
    L = 0.000987*h*(1 - (Zc1/Zc2)*sqrt(e_eff1/e_eff2))^2;

    L1 = (Lw1/(Lw1+Lw2))*L; % nH
    L2 = (Lw2/(Lw1+Lw2))*L; % nH
    c1 = sqrt(e_eff1)/Zc1;
    c2 = 1-(W2/W1);
    c3 = (e_eff1+0.3)/(e_eff1-0.258);
    c4 = ((W1/h)+0.264)/((W1/h)+0.8);
    C = 0.00137*h*c1*c2*c3*c4; % pF
    if cambio
        Laux = L1;
        L1 = L2;
        L2 = Laux;
    end
end

function [Y1_dB, Y3_dB] = scattering_fc(Zo, WL, WC, N, fc, h, er, l_L, l_C, f_c)
    % esta función sirve para calcular la atenuación exactamente en la fc,
    % pero utilizando una f_c un poco menor para que sean exactos -3dB
    c = 3e8; % velocidad de la luz
    [Z_0L, l_gL, e_effL] = Z_c(WL, fc, h, er);
    [Z_0C, l_gC, e_effC] = Z_c(WC, fc, h, er);
    z=[Z_0L Z_0C Z_0L Z_0C Z_0L Z_0C Z_0L Z_0C Z_0L];
    Z = z(1:N);
    e_e = [e_effL e_effC e_effL e_effC e_effL e_effC e_effL e_effC e_effL];
    e_eff = e_e(1:N);
    l = (l_L+l_C).*1e-3;
    
    betam =(2*pi*f_c*1e9/c).*sqrt(e_eff); % (2*pi*f*1e9/c).*sqrt(e_eff);

    A = cos(betam.*l); % cos(BL.*f./Fc2);
    B = 1i.*Z.*sin(betam.*l); % *sin(BL.*f./Fc2);
    C = 1i.*(1./Z).*sin(betam.*l); % sin(BL.*f./Fc2);
    D = cos(betam.*l); % cos(BL.*f./Fc2);
    %MM es una matriz de tamaño 2x2xN
    MM(1,1,:) = A;  MM(1,2,:) = B;
    MM(2,1,:) = C;  MM(2,2,:) = D;

    %Multiplicación Matricial
    T = eye(2);
    for Ind=1:N
        T = T*MM(:,:,Ind);
    end
    A0=T(1,1); B0=T(1,2); C0=T(2,1); D0=T(2,2); % Coeficientes A,B,C y D

    den = A0 + (B0/Zo) + (C0*Zo) + D0;
    S11 = (A0 + (B0/Zo) - (C0*Zo) - D0)/den;
    Y1_dB = 20*log10(abs(S11));

    S21 = 2/den; % (A0+(B0./Zo)+(C0.*Zo)+D0).\2;
    Y3_dB = 20*log10(abs(S21));
end

function [Y1_dB, Y3_dB, frec] = scattering(Zo, WL, WC, N, fc, h, er)
    c = 3e8; % velocidad de la luz

    freqmi=0*1e9;
    freqma=4*fc*1e9;
    freqstep=0.05*1e9;         
    nfreqs=((freqma-freqmi)/freqstep);
    freq = freqmi:freqstep:freqma; % freq(1:nfreqs)=freqmi:freqstep:freqma;

    W_0 = calc_W(Zo, er, h, fc); 
    [Z_0, l_g, e_eff0] = Z_c(W_0, fc, h, er);
    [Z_0L, l_gL, e_effL] = Z_c(WL, fc, h, er);
    [Z_0C, l_gC, e_effC] = Z_c(WC, fc, h, er);

    [l_L, l_C] = calc_largo(W_0, WL, WC, fc, N, h, er);

    z=[Z_0L Z_0C Z_0L Z_0C Z_0L Z_0C Z_0L Z_0C Z_0L];
    Z = z(1:N);
    l = (l_L+l_C).*1e-3;  
    e_e = [e_effL e_effC e_effL e_effC e_effL e_effC e_effL e_effC e_effL];
    e_eff = e_e(1:N);
    
    [L1, L2, C] = Step_Width(W_0, WL, Z_0, Z_0L, e_eff0, e_effL, h);
    [L1, L2, C] = Step_Width(WL, WC, Z_0L, Z_0C, e_effL, e_effC, h);

    n = 0;
    dif = -1;
    ant = 0;

    while (n<nfreqs)&&(dif < 0)
        n = n+1;
        frec(n) = freq(n);
        % ABCD según "2-D Electromagnetic Simulation of Passive Microstrip
        % Circuits" p 31 parece correcto
        %betac(K)=2*pi*f*sqrt(muz*epsz*epsrc);
        betam =(2*pi*freq(n)/c).*sqrt(e_eff); % (2*pi*f*1e9/c).*sqrt(e_eff);

        A=cos(betam.*l); % cos(BL.*f./Fc2);
        B=1i.*Z.*sin(betam.*l); % *sin(BL.*f./Fc2);
        C=1i.*(1./Z).*sin(betam.*l); % sin(BL.*f./Fc2);
        D=cos(betam.*l); % cos(BL.*f./Fc2);
        %MM es una matriz de tamaño 2x2xN
        MM(1,1,:) = A;
        MM(1,2,:) = B;
        MM(2,1,:) = C;
        MM(2,2,:) = D;

        %Multiplicación Matricial
        T = eye(2);
        for Ind=1:N
            T = T*MM(:,:,Ind);
        end
        A0=T(1,1); B0=T(1,2); C0=T(2,1); D0=T(2,2); % Coeficientes A,B,C y D
        
        den = A0 + (B0/Zo) + (C0*Zo) + D0;
        S11 = (A0 + (B0/Zo) - (C0*Zo) - D0)/den;
        Y1(n)=S11; 
        Y1_dB(n)=20*log10(abs(S11));
        S21 = 2/den; % (A0+(B0./Zo)+(C0.*Zo)+D0).\2;
        Y3(n)=S21; 
        Y3_dB(n)=20*log10(abs(S21));
        
        if (freq(n) > fc*1e9)
            dif = Y3_dB(n)-ant;
        end
        ant = Y3_dB(n);
    end
end

function [W] = calcW(Zc, e_r, h)
    A = (Zc/60)*((e_r+1)/2)^0.5 + ((e_r-1)/(e_r+1))*(0.23+(0.11/e_r));
    u = 8*exp(A)/(exp(2*A)-2);
    if (u>2)
        B = (60*pi^2)/(Zc*sqrt(e_r));
        u = (2/pi)*((B-1)-log(2*B-1)+((e_r-1)/(2*e_r))*(log(B-1)+0.39-(0.61/e_r)));
    end
    W = u*h;
end

function [W] = calc_W(Zc, er, h, fc)
    zc = Zc;
    W = calcW(Zc, er, h);
    Z = Z_c(W, fc, h, er);
    dif = Z-Zc;
    i = 0;
    while (abs(dif) > 0.0001)&&(i<4)
        zc = zc-dif;        
        W = calcW(zc, er, h);
        Z = Z_c(W, fc, h, er);
        dif = Z-Zc;
        i = i+1;
    end
end

function boton_graf_scattering(hObject, eventdata)
    global app
    Zo = str2double(get(app.e_Zo,'String'));
    WL = str2double(get(app.e_WL,'String'));
    WC = str2double(get(app.e_WC,'String'));
    N = str2double(get(app.e_n,'String'));
    fc = str2double(get(app.e_fc,'String'));
    h = str2double(get(app.e_h,'String'));
    er = str2double(get(app.e_eps_r,'String'));
    
    W0 = calc_W(Zo, er, h, fc); 
    [l_L, l_C, l_Zo, f_c] = calc_largo_fc(Zo, W0, WL, WC, fc, N, h, er)
        
    figure()
    [Y1_dB, Y2_dB, freq] = scattering(Zo, WL, WC, N, f_c, h, er);
    freqGHz = freq*1e-9;
    dfreq = (freqGHz(2)-freqGHz(1));
    d_Y2 = diff(Y2_dB)./dfreq;
    dY2 = [0 d_Y2];
    plot(freqGHz, Y1_dB, 'r', freqGHz, Y2_dB, 'k','LineWidth',1.5);
    grid on
%     axis([0 max(freq) -40 2])
    title('Scattering parameters','Fontsize',12)
    xlabel('Frequency, f (GHz)','Fontsize',12)
    ylabel('|S_{11}| and |S_{21}|, (dB)','Fontsize',12)
    dataSim(:,1) = freqGHz';
    dataSim(:,2) = Y1_dB';
    dataSim(:,3) = Y2_dB';
    save('datos_sim.mat', 'dataSim')
    
    figure()
    plot(freqGHz, -dY2, 'LineWidth',1.5);
    grid on
    xlabel('Frequency, f (GHz)','Fontsize',12)
    ylabel('L_A´, (dB/GHz)','Fontsize',12)
end

function boton_print_data(hObject, eventdata)
    global app
    Zo = str2double(get(app.e_Zo,'String'));
    WL = str2double(get(app.e_WL,'String'));
    WC = str2double(get(app.e_WC,'String'));
    N = str2double(get(app.e_n,'String'));
    fc = str2double(get(app.e_fc,'String'));
    h = str2double(get(app.e_h,'String'));
    er = str2double(get(app.e_eps_r,'String'));
    filterString = get(app.pm_TF,'String');
    filterValue = get(app.pm_TF,'Value');
        
    W_0 = calc_W(Zo, er, h, fc);
    [Z_0, l_g0] = Z_c(W_0, fc, h, er); %Impedancia_c(W_0, fc); 
    [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, WL, WC, fc, N, h, er);
    largo = 2*l_Z0 + sum(l_C+l_L);
    if imag(largo) ~= 0
        largo = 0;
        Min_der = 0;
        Min_S12 = 0;
    else
        [Y1_dB, Y2_dB, freq] = scattering(Zo, WL, WC, N, f_c, h, er);
        dfreq = (freq(2)-freq(1))*1e-9;
        d_Y2 = diff(Y2_dB)./dfreq;
        Min_der = min(d_Y2);
        Min_S12 = min(Y2_dB);
    end
    [lim_WL, lim_WC] = encuentra_limites(W_0, fc, N, h, er);

    g = param_g(N);
    L = zeros(1,N);
    C = L;
    Fc = fc*1e9;
    wc = 2*pi*Fc;
    for i = 1:N
        if (mod(i,2))  % i es impar?
            L(i) = (Z_0/wc)*g(i);
        else
            C(i) = (1/(Z_0*wc))*g(i);
        end
    end
        
    ancho = WC+2;
    dataLengths = {'   l_Zo' l_Z0 ' ' 0;
                   '   Wo' W_0 ' ' 0;
                   '   W_L' WL ' ' 0;
                   '   W_C' WC ' ' 0};
    for i = 1:N
       if (mod(i,2)) % (l_L(i) ~= 0)
           text = ['   l_L' num2str(i)];
           text2 = ['   L_' num2str(i)];
           dataLengths = [dataLengths; {text, l_L(i), text2, L(i)}];
       else
           text = ['   l_C' num2str(i)];
           text2 = ['   C_' num2str(i)];
           dataLengths = [dataLengths; {text, l_C(i), text2, C(i)}];   
       end
    end    
    
    a = 0;
    b = 0;
    if IsOctave()
        a = 35;
        b = 10;
    end
    d = figure('Name','Filter Data',...%'Color',[0.9400 0.9400 0.9400],...
               'MenuBar','none','NumberTitle','off',...
               'Position',[400 200 535-2*b 260+a]);
  
    % Tabla
    cnames = {'Parameter','Value'};
    cformat = {'char', 'numeric'};
    cedit = [true true];
    cwidth = {155 60};
    datosB = [filterString(filterValue), {0}];
    datos = {'Filter order, n', N;
             'Cutoff frequency, fc (GHz)', fc;
             'Source impedance, Zo (Ohms)' Z_0;
             'Dielectric constant, Eo' er;
             'Substrate height, h (mm)' h;
             'Limit of max W_L' lim_WL;
             'Limit of min W_C' lim_WC;
             'Total length (mm)', largo;
             'Total width (mm)', ancho;
             'max(L_A´), (dB/GHz)', -Min_der;
             'max(L_A), (dB)', -Min_S12};
    datos = [datosB; datos];
    ui_tabla = uitable('Position',[10,10,250-b,240+a],...
                 'FontName','Microsoft Sans Serif',...
                 'FontSize',8, 'Data',datos,...
                 'ColumnName',cnames,...
                 'ColumnWidth',cwidth,'ColumnEditable',cedit,...
                 'ColumnFormat',cformat);
     
     % Tabla 2
    cnames = {'Length','Val. (mm)', 'L, C', 'Val. (H, F)'};
    cformat = {'char', 'numeric','char', 'numeric'};
    cedit = [false false false false];
    cwidth = {50 65 43 65};
    ui_tabla = uitable('Position',[270-b,10,258-b,240+a],...
                 'FontName','Microsoft Sans Serif',...
                 'FontSize',8, 'Data',dataLengths,...
                 'ColumnName',cnames,...
                 'ColumnWidth',cwidth,'ColumnEditable',cedit,...
                 'ColumnFormat',cformat);
end

function [lim_WL, lim_WC] = encuentra_limites(W_0, fc, n, h, er)
    [Zo] = Z_c(W_0, fc, h, er);
    lim_WC = 20*W_0; % WC(i_WC);
    lim_WC2 = W_0; % WC(i_WC-1); % es el imaginario

    lim_WL = 0.01; % WL(i_WL-1);
    lim_WL2 = W_0; % WL(i_WL);  % es el imaginario
    
    precision = 1e-5; % 0.00001;    
    
    dif_WL = lim_WL2 - lim_WL;
    while (dif_WL > precision)
        WL_int = lim_WL + dif_WL/2;
        W_C = 20*W_0;
        [l_L, l_C] = calc_largo(W_0, WL_int, W_C, fc, n, h, er);  
        len = sum(l_C+l_L);
        if imag(len) ~= 0
            lim_WL2 = WL_int;
        else
            lim_WL = WL_int;
        end
        dif_WL = lim_WL2 - lim_WL;
    end
    
    dif_WC = lim_WC - lim_WC2;
    while (dif_WC > precision)
        WC_int = lim_WC2 + dif_WC/2;
        W_L = lim_WL/2; 
        [l_L, l_C] = calc_largo(W_0, W_L, WC_int, fc, n, h, er);
        len = sum(l_C+l_L);
        if imag(len) ~= 0
            lim_WC2 = WC_int;
        else
            lim_WC = WC_int;
        end
        dif_WC = lim_WC - lim_WC2;
    end
    
end

function boton_find_limits(~, ~)
    global app
    N = str2double(get(app.e_n,'String'));
    fc = str2double(get(app.e_fc,'String'));
    h = str2double(get(app.e_h,'String'));
    er = str2double(get(app.e_eps_r,'String'));      
    Zo = str2double(get(app.e_Zo,'String')); 
    W_0 = calc_W(Zo, er, h, fc);

    [lim_WL, lim_WC] = encuentra_limites(W_0, fc, N, h, er);
    
    set(app.e_WLmax,'String',num2str(lim_WL));
    set(app.e_WCmin,'String',num2str(lim_WC));
    %---------------------------------------------------------------------
    figure()
    hold on
    filtro = get(app.pm_TF,'Value');
    WLmin = str2double(get(app.e_WLmin,'String')); 
    WCmax = str2double(get(app.e_WCmax,'String')); 
    n = 3:9;
    for fil = 1:6        
        set(app.pm_TF,'Value',fil);
        for i = 1:length(n)
            W0(i) = W_0;
            [lim_WL(i, fil), lim_WC(i,fil)] = encuentra_limites(W_0, fc, n(i), h, er);
            if (fil > 2)
                WLmax = lim_WL(i, fil)*0.93; % *0.97;
                WCmin = lim_WC(i,fil)+0.0001;
               
                [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, WLmax, WCmin, fc, n(i), h, er);
                TLmax(i, fil-2) = 2*l_Z0 + sum(l_C+l_L);
                
                [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, WLmin, WCmax, fc, n(i), h, er);
                TLmin(i, fil-2) = 2*l_Z0 + sum(l_C+l_L);
            
                [Y1_dB, Y2_dB, freq] = scattering(Zo, WLmin, WCmax, n(i), f_c, h, er);
                dfreq = (freq(2)-freq(1))*1e-9;
                d_Y2 = diff(Y2_dB)./dfreq;
                LAp(1) = -min(d_Y2);   LA(1) = -min(Y2_dB); 
                
                [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, WLmax, WCmin+1, fc, n(i), h, er);
                [Y1_dB, Y2_dB, freq] = scattering(Zo, WLmax, WCmin+1, n(i), f_c, h, er);
                dfreq = (freq(2)-freq(1))*1e-9;
                d_Y2 = diff(Y2_dB)./dfreq;
                LAp(2) = -min(d_Y2);   LA(2) = -min(Y2_dB); 
                
                [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, WLmax, WCmax, fc, n(i), h, er);
                [Y1_dB, Y2_dB, freq] = scattering(Zo, WLmax, WCmax, n(i), f_c, h, er);
                dfreq = (freq(2)-freq(1))*1e-9;
                d_Y2 = diff(Y2_dB)./dfreq;
                LAp(3) = -min(d_Y2);   LA(3) = -min(Y2_dB);  
                
                [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, WLmin, WCmin+1, fc, n(i), h, er);
                [Y1_dB, Y2_dB, freq] = scattering(Zo, WLmin, WCmin+1, n(i), f_c, h, er);
                dfreq = (freq(2)-freq(1))*1e-9;
                d_Y2 = diff(Y2_dB)./dfreq;
                LAp(4) = -min(d_Y2);   LA(4) = -min(Y2_dB); 
                
                if (i<7)
                    LApmax(i, fil-2) = max(LAp);
                    LAmax(i, fil-2) = max(LA); 
                    LApmin(i, fil-2) = min(LAp);
                    LAmin(i, fil-2) = min(LA); 
                else
                    LAp(1)
                    LAp(3)
                    LApmax(i, fil-2) = max(LAp(1),LAp(3));
                    LAmax(i, fil-2) = max(LA(1),LA(3)); 
                    LApmin(i, fil-2) = min(LAp(1),LAp(3));
                    LAmin(i, fil-2) = min(LA(1),LA(3)); 
                end
            end
        end   
    end     
    
%     lim_WL
%     lim_WC    
    p1 = plot(n, lim_WL(:, 1), 'r-o', 'LineWidth', 1.2);
    p2 = plot(n, lim_WL(:, 2), 'g-*', 'LineWidth', 1.2);
    p3 = plot(n, lim_WL(:, 3), 'b-s', 'LineWidth', 1.2);
    p4 = plot(n, lim_WL(:, 4), 'k-.d', 'LineWidth', 1.2);
    p5 = plot(n, lim_WL(:, 5), 'm--^', 'LineWidth', 1.2);   
    p6 = plot(n, lim_WL(:, 6), 'b:p', 'LineWidth', 1.2);
    %plot(n, lim_WC, tipo)
    h1 = plot(n, lim_WC(:, 1), 'r-o', 'LineWidth', 1.2);
    h2 = plot(n, lim_WC(:, 2), 'g-*', 'LineWidth', 1.2);
    h3 = plot(n, lim_WC(:, 3), 'b-s', 'LineWidth', 1.2);
    h4 = plot(n, lim_WC(:, 4), 'k-.d', 'LineWidth', 1.2);
    h5 = plot(n, lim_WC(:, 5), 'm--^', 'LineWidth', 1.2);   
    h6 = plot(n, lim_WC(:, 6), 'b:p', 'LineWidth', 1.2);
    %W0
    pZo = plot(n,W0,'k--', 'LineWidth', 2);
    grid on
    %xlim([3 15])
    ymax = ceil(max(max(lim_WC)));
    ylim([0 ymax+1.5])
    xlabel('Filter order, {n}')
    ylabel('Limit width, (mm)')
    text(7, W_0+0.4, 'W_0')
    
    a = get(app.pm_TF,'String');
    %a(7) = {'Width W_0'};
    lgd = legend([p1 p2 p3 p4 p5 p6], a,... % [p1 p2 p3 p4 p5 p6 pZo],a,...
        'Location','northeast','NumColumns',2);
    
    set(app.pm_TF,'Value',filtro)
    
    %----------------------------------------------------------------------
    figure()
    hold on
    p1 = plot(n, TLmax(:, 1), 'r-o', 'LineWidth', 1.2);
    p2 = plot(n, TLmax(:, 2), 'b-*', 'LineWidth', 1.2);
    p3 = plot(n, TLmax(:, 3), 'g-s', 'LineWidth', 1.2);
    p4 = plot(n, TLmax(:, 4), 'k-.d', 'LineWidth', 1.2);

    h1 = plot(n, TLmin(:, 1), 'r-o', 'LineWidth', 1.2);
    h2 = plot(n, TLmin(:, 2), 'b-*', 'LineWidth', 1.2);
    h3 = plot(n, TLmin(:, 3), 'g-s', 'LineWidth', 1.2);
    h4 = plot(n, TLmin(:, 4), 'k-.d', 'LineWidth', 1.2);

    grid on
    xlabel('Filter order, {n}')
    ylabel('TL, (mm)')
    lgd = legend([p1 p2 p3 p4], a(3:6),... 
        'Location','northwest');

    %----------------------------------------------------------------------
    figure()
    hold on
    p1 = plot(n, LAmax(:, 1), 'r-o', 'LineWidth', 1.2);
    p2 = plot(n, LAmax(:, 2), 'b-*', 'LineWidth', 1.2);
    p3 = plot(n, LAmax(:, 3), 'g-s', 'LineWidth', 1.2);
    p4 = plot(n, LAmax(:, 4), 'k-.d', 'LineWidth', 1.2);

    h1 = plot(n, LAmin(:, 1), 'r-o', 'LineWidth', 1.2);
    h2 = plot(n, LAmin(:, 2), 'b-*', 'LineWidth', 1.2);
    h3 = plot(n, LAmin(:, 3), 'g-s', 'LineWidth', 1.2);
    h4 = plot(n, LAmin(:, 4), 'k-.d', 'LineWidth', 1.2);

    grid on
    xlabel('Filter order, {n}')
    ylabel('L_A, (dB)')
    lgd = legend([p1 p2 p3 p4], a(3:6),... 
        'Location','northwest');
    
    %----------------------------------------------------------------------
    figure()
    hold on
    p1 = plot(n, LApmax(:, 1), 'r-o', 'LineWidth', 1.2);
    p2 = plot(n, LApmax(:, 2), 'b-*', 'LineWidth', 1.2);
    p3 = plot(n, LApmax(:, 3), 'g-s', 'LineWidth', 1.2);
    p4 = plot(n, LApmax(:, 4), 'k-.d', 'LineWidth', 1.2);

    h1 = plot(n, LApmin(:, 1), 'r-o', 'LineWidth', 1.2);
    h2 = plot(n, LApmin(:, 2), 'b-*', 'LineWidth', 1.2);
    h3 = plot(n, LApmin(:, 3), 'g-s', 'LineWidth', 1.2);
    h4 = plot(n, LApmin(:, 4), 'k-.d', 'LineWidth', 1.2);

    grid on
    xlabel('Filter order, {n}')
    ylabel("L_A', (dB)")
    lgd = legend([p1 p2 p3 p4], a(3:6),... 
        'Location','northwest');
    
end

function boton_simular(hObject, eventdata)
    global app
    Zo = str2double(get(app.e_Zo,'String'));
    WL = str2double(get(app.e_WL,'String'));
    WC = str2double(get(app.e_WC,'String'));
    n = str2double(get(app.e_n,'String'));
    fc = str2double(get(app.e_fc,'String'));
    h = str2double(get(app.e_h,'String'));
    er = str2double(get(app.e_eps_r,'String'));   
    
    W_0 = calc_W(Zo, er, h, fc);
    app.W_0 = W_0;
    
    W_Lmax = str2double(get(app.e_WLmax,'String')); % 1;
    W_Lmin = str2double(get(app.e_WLmin,'String')); % 0.5;
    app.paso_WL = str2double(get(app.e_WLstep,'String')); % 0.01;
    WL = W_Lmin:app.paso_WL:W_Lmax

    W_Cmax = str2double(get(app.e_WCmax,'String')); % 50;
    W_Cmin = str2double(get(app.e_WCmin,'String')); % 4; % debe ser mayor que W_0
    app.paso_WC = str2double(get(app.e_WCstep,'String')); % 0.5;
    WC = W_Cmin:app.paso_WC:W_Cmax
    
    for i = 1:length(WL)
        for j = 1:length(WC)
            [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, WL(i), WC(j), fc, n, h, er);
            largo(i, j) = 2*l_Z0 + sum(l_C+l_L);
            if imag(largo(i, j)) ~= 0
                largo(i, j) = NaN;
                Min_der(i, j) = NaN;
                Min_S12(i, j) = NaN;
            else
                [Y1_dB, Y2_dB, freq] = scattering(Zo, WL(i), WC(j), n, f_c, h, er);
                dfreq = (freq(2)-freq(1))*1e-9;
                d_Y2 = diff(Y2_dB)./dfreq;
                Min_der(i, j) = min(d_Y2);
                Min_S12(i, j) = min(Y2_dB);
            end
            ancho(i,j) = WC(j)+2;        
        end       
        por100 = 100*i/length(WL);
        cadena = [num2str(round(por100)) ' %']
        set(app.t_porcentaje,'String',cadena);
        refresh
    end    
    
    app.WC = WC;
    app.WL = WL;
    app.Min_der = Min_der;
    app.Min_S12 = Min_S12;
    app.largo = largo;
    
    set(app.btn_bgtl,'Enable','on');
    set(app.btn_lconst,'Enable', 'on');
    
    maxLen = max(max(largo));
    minLen = min(min(largo));
    set(app.e_minLen,'String',num2str(minLen));
    set(app.e_maxLen,'String',num2str(maxLen));
end

function boton_graf_MinDer(~, ~)
    global app
    WC = app.WC;
    WL = app.WL;
    Min_der = app.Min_der;
    figure()

    mesh(WC, WL, -Min_der)
    zlabel("max(L_A')")
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
end

function boton_graf_MinS21(~, ~)
    global app
    WC = app.WC;
    WL = app.WL;
    Min_S12 = app.Min_S12;
    figure()

    mesh(WC, WL, -Min_S12)
    zlabel('max(L_A)')
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
end

function boton_graf_length(~, ~)
    global app
    WC = app.WC;
    WL = app.WL;
    largo = app.largo;
    figure()

    mesh(WC, WL, largo)
    zlabel('Total length, (mm)')
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
    
    boton_graf_MinS21();
    
    boton_graf_MinDer();
end

function boton_lconst(~, ~)
    global app
    WC = app.WC;
    WL = app.WL;
    Min_der = app.Min_der;
    Min_S12 = app.Min_S12;
    largo = app.largo;
    paso_WL = app.paso_WL;
    paso_WC = app.paso_WC;
    n = str2double(get(app.e_n,'String'));
    fc = str2double(get(app.e_fc,'String'));
    h = str2double(get(app.e_h,'String'));
    er = str2double(get(app.e_eps_r,'String'));  
    Zo = str2double(get(app.e_Zo,'String'));
    W_0 = app.W_0;
    [Z_0, l_g0] = Z_c(W_0, fc, h, er); 
   
    
    %% Encuentra superficie para Largo constante
    Largo_const = str2double(get(app.e_cLen,'String')); % 60;
    
    lCons = Largo_const*ones(size(largo));
    [ren col] = size(largo);
    for i = 1:ren
        for j = 1:col
            if (isnan(largo(i,j)))
                lCons(i,j) = NaN;
            end
        end
    end
            
    k = 1;
    for i_WL = 1:length(WL)-1
        for i_WC = 1:length(WC)-1
            if (largo(i_WL,i_WC)>Largo_const)&&(largo(i_WL,i_WC+1)<Largo_const)
                LM = largo(i_WL,i_WC);
                Lm = largo(i_WL,i_WC+1);
                difTot = LM-Lm;
                difCons = Largo_const-Lm;
                pC = difCons/difTot;
                xx(k) = WC(i_WC+1)-pC*paso_WC;
                yy(k) = WL(i_WL);
                 
                [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, yy(k), xx(k), fc, n, h, er);
                Largo = 2*l_Z0 + sum(l_C+l_L);
                eLargo = Largo_const - Largo;
                wc = xx(k);  wl = yy(k);
                z = 1;
                while (abs(eLargo) > 0.001)&&(z<10)
                    wc = wc - eLargo*0.4;
                    xx(k) = wc;
                    [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, wl, wc, fc, n, h, er);
                    Largo = 2*l_Z0 + sum(l_C+l_L);
                    eLargo = Largo_const - Largo;
                    z = z+1;
                end
                
                k = k+1;
            end
            if (largo(i_WL,i_WC)<Largo_const)&&(largo(i_WL+1,i_WC)>Largo_const)
                LM = largo(i_WL+1,i_WC);
                Lm = largo(i_WL,i_WC);
                difTot = LM-Lm;
                difCons = Largo_const-Lm;
                pC = difCons/difTot;
                xx(k) = WC(i_WC);
                yy(k) = WL(i_WL)+pC*paso_WL;
                [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, yy(k), xx(k), fc, n, h, er);
                
                Largo = 2*l_Z0 + sum(l_C+l_L);
                eLargo = Largo_const - Largo;
                wc = xx(k);  wl = yy(k);
                z = 1;
                while (abs(eLargo) > 0.001)&&(z<10)
                    wl = wl + eLargo*0.01;
                    yy(k) = wl;
                    [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, wl, wc, fc, n, h, er);
                    Largo = 2*l_Z0 + sum(l_C+l_L);
                    eLargo = Largo_const - Largo;
                    z = z+1;
                end
                
                k = k+1;
            end
        end
    end
    k = k-1;    

    for i = 1:k
        W_C = xx(i);
        W_L = yy(i); 
        [l_L, l_C, l_Z0, f_c] = calc_largo_fc(Zo, W_0, W_L, W_C, fc, n, h, er);

        [Y1_dB, Y2_dB, freq] = scattering(Z_0, W_L, W_C, n, f_c, h, er); 
        dfreq = (freq(2)-freq(1))*1e-9;
        d_Y2 = diff(Y2_dB)./dfreq;
        min_der(i) = min(d_Y2);
        min_S12(i) = min(Y2_dB);

        Largo(i) = 2*l_Z0 + sum(l_C+l_L);
        Ancho(i) = W_C+2;          
        Area(i) = Largo_const*Ancho(i);
    end

    figure()
    subplot(2,1,1)
    mesh(WC, WL, largo)%, 'FaceColor', 'flat')
    hold on
    mesh(WC, WL, lCons, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none')
    zlabel('length')
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
    subplot(2,1,2)
    plot(xx,yy,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
        
    % Total Width
    figure()
    subplot(2,2,1)    
    plot(xx,yy,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
    subplot(2,2,2)
    plot(xx, Ancho,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel('Total Width')
    subplot(2,2,3)
    plot(yy, Ancho,'LineWidth',1.5)
    grid on
    xlabel('W_L (mm)')
    ylabel('Total Width')
    subplot(2,2,4)
    plot3(xx,yy,Ancho,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
    zlabel('Total Width')
    min_Width = min(Ancho);
    p = find(Ancho == min_Width);
    W_C_opt_width = xx(p);
    W_L_opt_width = yy(p);
    
    % Max(L_a')
    figure()
    subplot(2,2,1) 
    plot(xx,yy,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
    subplot(2,2,2)
    plot(xx, -min_der,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel("max(L_A')")
    subplot(2,2,3)
    plot(yy, -min_der,'LineWidth',1.5)
    grid on
    xlabel('W_L (mm)')
    ylabel("max(L_A')")
    subplot(2,2,4)
    plot3(xx,yy,-min_der,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
    zlabel("max(L_A')")
    max_LAp = max(-min_der);
    p = find(-min_der == max_LAp);
    W_C_opt_L_Ap = xx(p);
    W_L_opt_L_Ap = yy(p);
    
    % Max L_A
    figure()
    subplot(2,2,1) 
    plot(xx,yy,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
    subplot(2,2,2)
    plot(xx, -min_S12,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel("max(L_A), dB")
    subplot(2,2,3)
    plot(yy, -min_S12,'LineWidth',1.5)
    grid on
    xlabel('W_L (mm)')
    ylabel("max(L_A), dB")
    subplot(2,2,4)
    plot3(xx,yy,-min_S12,'LineWidth',1.5)
    grid on
    xlabel('W_C (mm)')
    ylabel('W_L (mm)')
    zlabel("max(L_A), dB")    
    max_LA = max(-min_S12);
    p = find(-min_S12 == max_LA);
    W_C_opt = xx(p);
    W_L_opt = yy(p);
    
    
    d = figure('Name','Optimal Values',...%'Color',[0.9400 0.9400 0.9400],...
               'MenuBar','none','NumberTitle','off',...
               'Position',[400 200 340 110]);
  
    % Tabla
    cnames = {'Parameter','Value','W_C','W_L'};
    cformat = {'char', 'numeric', 'numeric', 'numeric'};
    cedit = [true true true true];
    cwidth = {100 60 60 60};
    datos = {'Min(Total Width)', min_Width, W_C_opt_width, W_L_opt_width;
             'Max(L_A)', max_LA, W_C_opt, W_L_opt;
             'Max(L_A´)', max_LAp, W_C_opt_L_Ap, W_L_opt_L_Ap};
    ui_tabla = uitable('Position',[10,10,320,90],...
                 'FontName','Microsoft Sans Serif',...
                 'FontSize',8, 'Data',datos,...
                 'ColumnName',cnames,...
                 'ColumnWidth',cwidth,'ColumnEditable',cedit,...
                 'ColumnFormat',cformat);             
end

function boton_lim_Zo(~,~)
    global app
    n = str2double(get(app.e_n,'String'));
    fc = str2double(get(app.e_fc,'String'));
    h = str2double(get(app.e_h,'String'));

    Zo_ini = 25;
    Zo_fin = 50;
    Zo_paso = 2.5;
    Zo = Zo_ini:Zo_paso:Zo_fin;
    er_ini = 4;
    er_paso = 2;
    er_fin = 12;
    er = er_ini:er_paso:er_fin;
    for i = 1:length(Zo)
        for j = 1:length(er)
            W_0(i,j) = calc_W(Zo(i), er(j), h, fc);
            [lim_WL(i, j), lim_WC(i, j)] = encuentra_limites(W_0(i, j), fc, n, h, er(j));
        end
        por100 = 100*i/length(Zo);
        disp([num2str(por100) ' %']);
    end
    
    figure()
    plot(Zo, W_0,'LineWidth',1.5)
    grid on
    xlabel('Z_0, (Ohms)','FontSize',12)  % ,'Interpreter','latex'
    ylabel('W_0, (mm)','FontSize',12) % ,'Interpreter','latex'
    for i = 1:length(er)
        texto = ['\epsilon_r = ' num2str(er(i))];
        text(Zo(5), W_0(5,i),texto, 'BackgroundColor','w')
    end
    
    figure()  
    subplot(2,1,1)
    plot(Zo, lim_WC,'LineWidth',1.5)
    grid on
    ylabel('Limit W_C, (mm)','FontSize',12) % ,'Interpreter','latex'
    for i = 1:length(er)
        texto = ['\epsilon_r = ' num2str(er(i))];
        text(Zo(i+1)+1, lim_WC(i+1,i),texto, 'HorizontalAlignment','l')
    end
        
    subplot(2,1,2)
    plot(Zo, lim_WL,'LineWidth',1.5)  
    grid on
    xlabel('Z_0, (Ohms)','FontSize',12)  % ,'Interpreter','latex'
    ylabel('Limit W_L, (mm)','FontSize',12) % ,'Interpreter','latex'
    for i = 1:length(er)
        texto = ['\epsilon_r = ' num2str(er(i))];
        text(Zo(i+1)+1, lim_WL(i+1,i),texto, 'HorizontalAlignment','l')
    end
end

function boton_length_fc(~,~)
    global app
    WC = str2double(get(app.e_WC,'String'));
    WL = str2double(get(app.e_WL,'String'));
    n = str2double(get(app.e_n,'String'));
    Zo = str2double(get(app.e_Zo,'String'));
    h = str2double(get(app.e_h,'String'));
    er = str2double(get(app.e_eps_r,'String')); 
    
    fc_ini = 0.5;
    fc_paso = 0.01;
    fc_fin = 5;
    fc = fc_ini:fc_paso:fc_fin;
    for i = 1:length(fc)
        W_0 = calc_W(Zo, er, h, fc(i));
        Z_0 = Z_c(W_0, fc(i), h, er); 

        [l_L, l_C, l_Z0] = calc_largo(W_0, WL, WC, fc(i), n, h, er);  
        largo(i) = 2*l_Z0 + sum(l_C+l_L);
        [Y1_dB, Y2_dB, freq] = scattering(Z_0, WL, WC, n, fc(i), h, er); 
        dfreq = (freq(2)-freq(1))*1e-9;
        d_Y2 = diff(Y2_dB)./dfreq;
    end
    figure()
    plot(fc, largo,'LineWidth',1.5)
    grid on
    xlabel('f_c, (GHz)','FontSize',12)  % ,'Interpreter','latex'
    ylabel('Total length, (mm)','FontSize',12) % ,'Interpreter','latex'
    
end

function isOctave = IsOctave()
% Returns true if this code is being executed by Octave.
% Returns false if this code is being executed by MATLAB, or any other MATLAB
% variant.
%
%    usage: isOctave = IsOctave()
    persistent octaveVersionIsBuiltIn;
    if (isempty(octaveVersionIsBuiltIn))
        octaveVersionIsBuiltIn = (exist('OCTAVE_VERSION', 'builtin') == 5);
        % exist returns 5 to indicate a built-in function.
    end
    isOctave = octaveVersionIsBuiltIn;
    % If OCTAVE_VERSION is a built-in function, then we must be in Octave.
    % Since the result cannot change between function calls, it is cached in a
    % persistent variable. isOctave cannot be a persistent variable, because it
    % is the return value of the function, so instead the persistent result must
    % be cached in a separate variable.
end