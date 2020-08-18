clc
close all
clearvars

% NEZNÁMÁ POLOHA STØEDU ROTACE
RX = 1;
RY = 1;

% NEZNÁMÉ ROZMÌRY DÍRY
l = 1;
d = 0.5;

% NEZNÁMÝ ÚHEL VYCHÝLENÍ
alpha = 20; % DEG
alpha = alpha * pi/180; %to RAD

% výchozí poloha drátu
wireX = [-2*d/(2*tan(abs(alpha))) 2*d/(2*tan(abs(alpha)))];
wireY = [wireX(1)*sin(alpha) wireX(2)*sin(alpha)];

% poloha horní hrany
upperX = [-l/2 l/2];
upperY = [d/2 d/2];

% poloha dolní hrany
lowerX = [-l/2 l/2];
lowerY = [-d/2 -d/2];

% vykreslení obrázku
figure('Name','Horní pohled');
graf = plot(RX, RY, 'r+', upperX, upperY, 'k', lowerX, lowerY, 'k', wireX, wireY, 'b');
xlabel('Osa X');
ylabel('Osa Y');
axis equal;
pause(0.5);

%% Lineární pohyby pro zjištìní úhlu alfa
% rozmìr a
a = d/2 - l/2*sin(abs(alpha));

fprintf('Pohyb ve smìru osy Y+ \n');
for ii = 1:10
    wireY = wireY - a/10; 
    graf(4).YData = wireY;
    pause(0.1)
end

fprintf('Pohyb ve smìru osy Y- \n');
for ii = 1:10
    wireY = wireY + a/5; 
    graf(4).YData = wireY;
    pause(0.1)
end

fprintf('Rozmìr a = %d m \n',2*a);

% rozmìr b
b = sign(alpha)*(d/tan(abs(alpha))) - l;

fprintf('Pohyb ve smìru osy X- (X+) \n');
for ii = 1:10
    wireX = wireX + b/10; 
    graf(4).XData = wireX;
    pause(0.1)
end

%% Rotace o zjištìný úhel alfa
% Výpoèet úhlu alfa
alphaKnown = sign(alpha) * atan((2*a)/(b)); % RAD
alphaKnown = alphaKnown *180/pi;% to DEG

fprintf('Úhel alfa je %d° \n', alphaKnown);
fprintf('Rotace kolem osy Z \n');

parallel = 0;
pootoceno = 0;

while ~parallel
    % Souøadnice bodu doteku
    % dle http://www.ambrsoft.com/TrigoCalc/Circles2/circlrLine_.htm
    u1 = [wireX(2)-wireX(1) wireY(2)-wireY(1)];
    n1 = [-u1(2) u1(1)];
    c1 = -n1(1)*wireX(1) -n1(2)*wireY(1);

    m = (-n1(1))/n1(2);
    d1 = -c1/n1(2);
    
    % vzdálenost hrany díry od drátu
    if alpha > 0 && RX >= 0
        rohX = upperX(2);
        rohY = upperY(2);
    elseif alpha > 0 && RX < 0
        rohX = lowerX(1);
        rohY = lowerY(1);
    elseif alpha < 0 && RX >= 0
        rohX = lowerX(2);
        rohY = lowerY(2);
    elseif alpha < 0 && RX < 0
        rohX = upperX(1);
        rohY = upperY(1);
    end

    r = sqrt((RX-rohX)^2 + (RY-rohY)^2);

    gama = r*r*(1+m*m) - (RY-m*RX-d1)^2;

    x1 = (RX + RY*m - d1*m + sqrt(gama))/(1+m*m);
    y1 = (d1 + RX*m + RY*m*m + m*sqrt(gama))/(1+m*m);

    x2 = (RX + RY*m - d1*m - sqrt(gama))/(1+m*m);
    y2 = (d1 + RX*m + RY*m*m - m*sqrt(gama))/(1+m*m);

    
    if imag(x2) == 0 % dojde k nárazu
        % Úhel doteku
        % dle https://www.mathworks.com/matlabcentral/answers/331017-calculating-angle-between-three-points
        P0 = [RX, RY];
        P1 = [rohX, rohY];
        if RX >= 0
            P2 = [x2, y2];
        else
            P2 = [x1, y1];
        end
        n1 = (P2 - P0) / norm(P2 - P0);  % Normalized vectors
        n2 = (P1 - P0) / norm(P1 - P0);

        angle3 = atan2(norm(det([n2; n1])), dot(n1, n2));
    else % nedojde k nárazu
        angle3 = inf;
    end

    if angle3 < abs(alpha - pootoceno)
        uhelPootoceni = - sign(alpha) * angle3;
        pootoceno = pootoceno - uhelPootoceni;
    else
        uhelPootoceni = abs(alpha - pootoceno);
        uhelPootoceni = - sign(alpha)*abs(alpha - pootoceno);
        parallel = 1;
    end

    uhelPootoceni = uhelPootoceni/10;
    % Rotace
    for jj = 1:10
        % Rotace horní úseèky
        for ii = 1:2
        vzdalenost = sqrt((upperX(ii)-RX)^2+(upperY(ii)-RY)^2);
        Beta = atan(abs(upperX(ii)-RX)/abs(upperY(ii)-RY));

            if RY <= upperY(ii) && RX >= upperX(ii)
                upperX(ii) = RX + vzdalenost * sin(-Beta + uhelPootoceni);
                upperY(ii) = RY + vzdalenost * cos(-Beta + uhelPootoceni);
            elseif RY <= upperY(ii) && RX <= upperX(ii)
                upperX(ii) = RX - vzdalenost * sin(-Beta - uhelPootoceni);
                upperY(ii) = RY + vzdalenost * cos(-Beta - uhelPootoceni);
            elseif RY >= upperY(ii) && RX <= upperX(ii)
                upperX(ii) = RX + vzdalenost * sin(Beta - uhelPootoceni);
                upperY(ii) = RY - vzdalenost * cos(Beta - uhelPootoceni);
            elseif RY >= upperY(ii) && RX >= upperX(ii)
                upperX(ii) = RX - vzdalenost * sin(Beta + uhelPootoceni);
                upperY(ii) = RY - vzdalenost * cos(Beta + uhelPootoceni);
            end
        end

        graf(2).XData = upperX;
        graf(2).YData = upperY;

        % Rotace dolní úseèky
        for ii = 1:2
        vzdalenost = sqrt((lowerX(ii)-RX)^2+(lowerY(ii)-RY)^2);
        Beta = atan(abs(lowerX(ii)-RX)/abs(lowerY(ii)-RY));

            if RY <= lowerY(ii) && RX >= lowerX(ii)
                lowerX(ii) = RX + vzdalenost * sin(-Beta + uhelPootoceni);
                lowerY(ii) = RY + vzdalenost * cos(-Beta + uhelPootoceni);
            elseif RY <= lowerY(ii) && RX <= lowerX(ii)
                lowerX(ii) = RX - vzdalenost * sin(-Beta - uhelPootoceni);
                lowerY(ii) = RY + vzdalenost * cos(-Beta - uhelPootoceni);
            elseif RY >= lowerY(ii) && RX <= lowerX(ii)
                lowerX(ii) = RX + vzdalenost * sin(Beta - uhelPootoceni);
                lowerY(ii) = RY - vzdalenost * cos(Beta - uhelPootoceni);
            elseif RY >= lowerY(ii) && RX >= lowerX(ii)
                lowerX(ii) = RX - vzdalenost * sin(Beta + uhelPootoceni);
                lowerY(ii) = RY - vzdalenost * cos(Beta + uhelPootoceni);
            end

        end

        graf(3).XData = lowerX;
        graf(3).YData = lowerY;
        pause(0.1)
    end
    
    if parallel == 1 %srovnání chyby vizualizace
        drat = atan((wireY(2)-wireY(1))/(wireX(2)-wireX(1)));
        dira = atan((upperY(2)-upperY(1))/(upperX(2)-upperX(1)));

        uhelPootoceni = dira - drat;
        %Rotace
        for jj = 1:1
            % Rotace horní úseèky
            for ii = 1:2
            vzdalenost = sqrt((upperX(ii)-RX)^2+(upperY(ii)-RY)^2);
            Beta = atan(abs(upperX(ii)-RX)/abs(upperY(ii)-RY));

                if RY <= upperY(ii) && RX >= upperX(ii)
                    upperX(ii) = RX + vzdalenost * sin(-Beta + uhelPootoceni);
                    upperY(ii) = RY + vzdalenost * cos(-Beta + uhelPootoceni);
                elseif RY <= upperY(ii) && RX <= upperX(ii)
                    upperX(ii) = RX - vzdalenost * sin(-Beta - uhelPootoceni);
                    upperY(ii) = RY + vzdalenost * cos(-Beta - uhelPootoceni);
                elseif RY >= upperY(ii) && RX <= upperX(ii)
                    upperX(ii) = RX + vzdalenost * sin(Beta - uhelPootoceni);
                    upperY(ii) = RY - vzdalenost * cos(Beta - uhelPootoceni);
                elseif RY >= upperY(ii) && RX >= upperX(ii)
                    upperX(ii) = RX - vzdalenost * sin(Beta + uhelPootoceni);
                    upperY(ii) = RY - vzdalenost * cos(Beta + uhelPootoceni);
                end
            end

            graf(2).XData = upperX;
            graf(2).YData = upperY;

            % Rotace dolní úseèky
            for ii = 1:2
            vzdalenost = sqrt((lowerX(ii)-RX)^2+(lowerY(ii)-RY)^2);
            Beta = atan(abs(lowerX(ii)-RX)/abs(lowerY(ii)-RY));

                if RY <= lowerY(ii) && RX >= lowerX(ii)
                    lowerX(ii) = RX + vzdalenost * sin(-Beta + uhelPootoceni);
                    lowerY(ii) = RY + vzdalenost * cos(-Beta + uhelPootoceni);
                elseif RY <= lowerY(ii) && RX <= lowerX(ii)
                    lowerX(ii) = RX - vzdalenost * sin(-Beta - uhelPootoceni);
                    lowerY(ii) = RY + vzdalenost * cos(-Beta - uhelPootoceni);
                elseif RY >= lowerY(ii) && RX <= lowerX(ii)
                    lowerX(ii) = RX + vzdalenost * sin(Beta - uhelPootoceni);
                    lowerY(ii) = RY - vzdalenost * cos(Beta - uhelPootoceni);
                elseif RY >= lowerY(ii) && RX >= lowerX(ii)
                    lowerX(ii) = RX - vzdalenost * sin(Beta + uhelPootoceni);
                    lowerY(ii) = RY - vzdalenost * cos(Beta + uhelPootoceni);
                end
            end

            graf(3).XData = lowerX;
            graf(3).YData = lowerY;
            pause(0.1)
        end
    end
    
    % výpoèet pohybu na doraz Y
    u1 = [wireX(2)-wireX(1) wireY(2)-wireY(1)];
    n1 = [-u1(2) u1(1)];
    c1 = -n1(1)*wireX(1) -n1(2)*wireY(1);

    if alpha > 0 && RX >= 0
        yTemp = (-n1(1)*lowerX(1) - c1)/n1(2);
        deltaTemp = lowerY(1) - yTemp;
    elseif alpha > 0 && RX < 0
        yTemp = (-n1(1)*upperX(2) - c1)/n1(2);
        deltaTemp = upperY(2) - yTemp;
    elseif alpha < 0 && RX >= 0
        yTemp = (-n1(1)*upperX(1) - c1)/n1(2);
        deltaTemp = upperY(1) - yTemp;
    elseif alpha < 0 && RX < 0
        yTemp = (-n1(1)*lowerX(2) - c1)/n1(2);
        deltaTemp = lowerY(2) - yTemp;
    end

    fprintf('Pohyb ve smìru osy Y \n');
    for ii = 1:10
        wireY = wireY + deltaTemp/10; 
        graf(4).YData = wireY;
        pause(0.1)
    end
end

%% Ustanovení ve støedové poloze na základì poloh dotekù
% dle
% https://www.e-matematika.cz/stredni-skoly/jak-urcit-obecnou-rovnici-primky-urcene-dvema-body.php
% a https://matematika.cz/vzdalenost-bod-primka

% vzdálenost horní hrany díry od drátu
u1 = [wireX(2)-wireX(1) wireY(2)-wireY(1)];
n1 = [-u1(2) u1(1)];
c1 = -n1(1)*wireX(1) -n1(2)*wireY(1);
pohyb1 = abs(n1(1)*lowerX(1) + n1(2)*lowerY(1) + c1)/(sqrt(n1(1)*n1(1) + n1(2)*n1(2)));

fprintf('Pohyb ve smìru osy Y+ \n');
for ii = 1:10
    wireY = wireY - pohyb1/10; 
    graf(4).YData = wireY;
    pause(0.1)
end

vPrev = wireY(1);
c1 = -n1(1)*wireX(1) -n1(2)*wireY(1);
pohyb2 = abs(n1(1)*upperX(1) + n1(2)*upperY(1) + c1)/(sqrt(n1(1)*n1(1) + n1(2)*n1(2)));

fprintf('Pohyb ve smìru osy Y- \n');
for ii = 1:10
    wireY = wireY + pohyb2/10; 
    graf(4).YData = wireY;
    pause(0.1)
end

pohyb2 = wireY(1) - vPrev;

fprintf('Pohyb ve smìru osy Y+ \n');
for ii = 1:5
    wireY = wireY - pohyb2/10; 
    graf(4).YData = wireY;
    pause(0.1)
end
pause(2)

%% Èelní pohled

% NEZNÁMÝ ÚHEL BETA
Beta = 10; %DEG
Beta = Beta * pi/180; %to RAD
wireZ = [-l/2*sin(Beta) l/2*sin(Beta)];

% vygenerování hrany díry
ang=0:0.01:2*pi; 
xp=d/2*cos(ang);
yp=d/2*sin(ang);

figure('Name','Èelní pohled');
graf2 = plot(xp, yp, 'k', [0 0], wireZ, '-o');
xlabel('Osa Y');
ylabel('Osa Z');
axis equal;

pohybZ = d/2 - abs(wireZ(2));

fprintf('Pohyb ve smìru osy Z+ \n');
for ii = 1:10
    wireZ = wireZ - pohybZ/10; 
    graf2(2).YData = wireZ;
    pause(0.1)
end

fprintf('Pohyb ve smìru osy Z- \n');
for ii = 1:10
    wireZ = wireZ + pohybZ/5; 
    graf2(2).YData = wireZ;
    pause(0.1)
end

fprintf('Pohyb ve smìru osy Y- (Y+) \n');
for ii = 1:10
    wireZ = wireZ - pohybZ/5; 
    graf2(2).YData = wireZ;
    pause(0.1)
end

fprintf('Návrat do výchozí polohy \n');
for ii = 1:10
    wireZ = wireZ + pohybZ/10; 
    graf2(2).YData = wireZ;
    pause(0.1)
end

fprintf('Úhel beta je %d° \n', Beta *180/pi);

vyrovnani = wireZ(2);

fprintf('Rotace kolem osy X \n');
for ii = 1:10
    wireZ(1) = wireZ(1) + vyrovnani/10;
    wireZ(2) = wireZ(2) - vyrovnani/10;
    graf2(2).YData = wireZ;
    pause(0.1)
end

fprintf('Pohyb ve smìru osy Z+ \n');
for ii = 1:10
    wireZ = wireZ - (d/2)/10; 
    graf2(2).YData = wireZ;
    pause(0.1)
end

fprintf('Pohyb ve smìru osy Z- \n');
for ii = 1:20
    wireZ = wireZ + (d/2)/10; 
    graf2(2).YData = wireZ;
    pause(0.05)
end

fprintf('Pohyb ve smìru osy Z+ \n');
for ii = 1:10
    wireZ = wireZ - (d/2)/10; 
    graf2(2).YData = wireZ;
    pause(0.1)
end