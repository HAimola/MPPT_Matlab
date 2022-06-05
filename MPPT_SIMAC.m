[Pow, Vs] = eval_Painel_Solar(298, 5, true, false);
[coef, S] = interpolar(Pow, Vs, 7, 'ponta_lin');
polyval(coef, 45.250536015760602) 

deriv = polyder(coef);
syms X
eq = poly2sym(deriv, X);
raizes = vpasolve(eq == 0);


root = inf;
resposta = 0;

for i = 1:length(raizes)
    if isreal(raizes(i))
      if polyval(deriv, double(raizes(i))) < root
         resposta = double(raizes(i));
      end
    end
end

sprintf("A resposta encontrada é %.15f, com um erro de %.15f",resposta, ...
       polyval(deriv, resposta))

function index = estrategia_de_interpolacao(estrategia, Vs, n)

    if isequal(estrategia, 'aleatorio')
        index = zeros(1, n+1);
        for m = 1:length(index)
            a = randi([1, 50]);
            
            if(~ismember(a, index)) 
                index(m) = a;
            else
                m = m - 1;
            end
        end
                   
        
    elseif isequal(estrategia, 'linear')
        index = int32(linspace(1, length(Vs)-1, n+1));
    
        
    elseif isequal(estrategia, 'ponta_lin')
        index = int32(linspace(length(Vs)/2, length(Vs), n+1));
    end
end

function [coef, S] = interpolar(Pow, Vs, n, estrategia)
    index = estrategia_de_interpolacao(estrategia, Vs, n);
    x = zeros(1, n+1);
    y = zeros(1, n+1);
    
    for i = 1:length(index)
        y(i) = Pow(index(i));
        x(i) = Vs(index(i));
    end

    [coef, S] = polyfit(x, y, n);
end


function [Currs, Vs] =  eval_Painel_Solar(Tin, resolution, PV, make_plot)

    % Constantes
    digits(15);
    k = 1.380649e-23; 
    q = 1.602176634e-19;
    Eg = 1.12; % Bandgap do Si

    % Constantes do Painel
    Ns = 12;
    Isc = 13.79;
    Voc = 49.34/Ns;

    N = 1.4; % Ideality factor
    Kisc = 0.048; % %/oC
    Rs = 0.08; % Resistencia em série

    syms V I G T

    Iph = (Isc - Kisc * (T - 298))*G/1000;
    Vt = k*T/q;
    Vt_ref = subs(Vt, T, 298);

    Io_ref = Isc/(exp(Voc/(N*Vt_ref))-1);
    Io = Io_ref * (T/298)^(3/N) * exp(-(Eg * q/(N*k)) * (1/T - 1/298));
    Ic = I - (Iph - Io*(exp((V+I.*Rs)/(Vt*N*Ns)) - 1))  == 0;

    irradiance = 250:250:1000;
    Vs = 0:1/resolution:Voc*Ns;
    
    
    for m = 1:length(irradiance)
        
        eq = subs(Ic, T, Tin);  
        eq = subs(eq, V, Vs);
        eq = subs(eq, G, irradiance(m));
        Currs = resolver_corrente(eq, Vs);
        
        if(PV)
            Currs = Currs .* Vs;
        end
    
     if(make_plot)
        hold on
        plot(Vs, Currs, 'o-', 'MarkerSize', 5, ...
            'MarkerIndices', 1:resolution*2:length(Vs), 'LineWidth', 1.3);
        
        txt = ['G = ' num2str(irradiance(m)) ' W/m^2'];
        
        if(~PV)
            text(Voc*Ns*0.5,Currs(1)+1, txt)
            xlim([0, Voc*Ns*1.05]);
            ylim([0, Currs(1)+3]);
        else
            text(38, m*125, txt, 'FontSize', 9)
            xlim([0, Voc*Ns*1.05]);
            ylim([0, 650]);
        end
        grid on
        grid minor

     end
    end

    if (make_plot)
        if(PV)
            tl = ['Curva P-V. Isc = ' num2str(Isc) '; Voc = ' ...
                                    num2str(Voc*Ns) '; N = ' num2str(N)];     
            ylabel("Potência de saída [W]"); 
        else
            tl = ['Curva I-V. Isc = ' num2str(Isc) '; Voc = ' ...
                                   num2str(Voc*Ns) '; N = ' num2str(N)];
            ylabel("Corrente de saída [A]");
        end
        title(tl);
        xlabel("Tensão sobre o módulo PV [V]");
    end  
end


function Currs = resolver_corrente(eq, Vs)
    syms I
    Currs = zeros(1, length(Vs));
    for c = 1:length(Vs)
        Currs(c) = vpasolve(eq(1, c), I);
    end 
end
