%Play Song
%Passion(flamma<->felix)
%Optimam(non-negans<->natalis)
%Rarus(determinans<->dies)
%Talentum(thensaurus<->tu)
%Armo(Alice<->mi)
%12+8, 19+17

%%% Calculate Third Octave Bands (base 2) in Matlab
%fcentre  = 10^3 * (2 .^ ([-18:13]/3))
%fd = 2^(1/6);
%fupper = fcentre * fd
%filedirectory = 'D:\PhD in Oxford,Ethz,KI,others\OxfordPhD\digital sound music(C++,matlab)\Frequency_Modulation_for_Digital_Synthesis_in_MATLAB\HornsE04.wav'

%[y, sr] = audioread('HornsE04Mono.wav'); 
display('Merry Christmas!')
%[y, sr] = audioread(fullfile(filedirectory));
%132300 = 44100*3
%sound(y, sr)
%pulsesep(1:pi/100:100*pi,10000);
%A frequency ratio expressed in octaves is the base-2 logarithm (binary logarithm) of the ratio:

%\log _{2}{\frac {f_{2}}{f_{1}}}
%An amplifier or filter may be stated to have a frequency response of �6 dB per octave over a particular frequency range, 
%which signifies that the power gain changes by �6 decibels (a factor of 4 in power), 
%when the frequency changes by a factor of 2.
%This slope, or more precisely 10 log10(4) ? 6.0206 decibels per octave, 
%corresponds to an amplitude gain proportional to frequency, which is equivalent to �20 dB per decade 
%(factor of 10 amplitude gain change for a factor of 10 frequency change). This would be a first-order filter.

%Ex: The distance between the frequencies 20 Hz and 40 Hz is 1 octave. An amplitude of 52 dB at 4 kHz decreases as frequency increases at ?2 dB/oct. What is the amplitude at 13 kHz?
%{\displaystyle {\text{number of octaves}}=\log _{2}\left({\frac {13}{4}}\right)=1.7}
%{\displaystyle {\text{Mag}}_{13{\text{ kHz}}}=52{\text{ dB}}+(1.7{\text{ oct}}\times -2{\text{ dB/oct}})=48.6{\text{ dB}}.\,}{\displaystyle {\text{Mag}}_{13{\text{ kHz}}}=52{\text{ dB}}+(1.7{\text{ oct}}\times -2{\text{ dB/oct}})=48.6{\text{ dB}}.\,}


%note = ['5','5','6','-','5','-','8','-';'7','-','-','-','-','-','-','-';'5','5','6','-','5','-','9','-';'8','-','-','-','-','-','-','-';];
startnoteoffset = 0;%middleC: +3
isminor = 0;
sr = 44100;%per second
speed = 60;
ns =2;
s = 8;%
ds = 1;
majors=[0 2 4 5 7 9 11 12];%startoffset+3 : middle major C
minors=[0 2 3 5 7 8 10 12];
if(isminor == 1)
    scale = majors;
else
    scale = minors;
end
%the statement above is equivalent to 
sr = sr/ds;
t = 1:sr;
t = t/sr;
%unit = sr/ns*speed/60*s/8;
%fundamental 220Hz, Harmonics: 440Hz,  660Hz...
%Octave Example on a guitar:
%from 0 (2: 110)
%1 octave(1/2: 55), 2 octaves(4:440)
%1. E (full length; 66 cm)
%2. E' (1/2 of length; one octave up; 164 Hz; 33 cm; touching the string in the middle suppress the fundamental or first harmonic so you hear this pitch)
%3. B' (1/3 of length; an octave and a fifth up; 246 Hz; 22 cm; touching the string here suppresses the second harmonic so you hear this pitch)
%4. E'' (1/4 of length; two octaves up; 328 Hz; 16.5 cm; touching the string here suppresses the third harmonic so you hear this pitch)
%5. G''' (1/5 of length; two octaves and a third up; 410 Hz; 13.2 cm; touching the string here suppresses the fourth harmonic so you hear this pitch; in fact a tempered G''' is 415.305 Hz, so the tuner will register this note as "sharp" if you touch the string right above the fret)

startnote = 220*(2^(startnoteoffset+3));
scale = startnote * scale;
%note = ['8', '7', '8', '6', '3', '7', '5', '3', '-'; '8', '8', '7', '6', '8', '5', '7', '5'];
%note =['6', '-', '7', '5', '4', '3', '2', '1', '7', '1', '7', '-', '6', '5', '-'];
PLOT = 1;
%TTLS = MakeScale(1,startnoteoffset, isminor, sr, ds, s, '*');
line = 1;
Note = zeros(8,8);
Note(1,:) = [5,5,6,6,5,5,9,9];
Note(2,:) = [7,7,7,7,7,7,7,7];
Note(3,:) = [5,5,6,6,5,5,10,10];
Note(4,:) = [9,9,9,9,9,9,9,9];
Note(5,:) = [5,5,13,13,11,11,9,9];
Note(6,:) = [7,7,6,6,6,6,6,6];
Note(7,:) = [12,12,11,11,9,9,10,10];
Note(8,:) = [9,9,9,9,9,9,9,9];

Outarray = zeros(size(Note, 1), sr*s);
%outarray = zeros(1,sr*s);
outarray = zeros(s,sr);
outarraymid = zeros(s,sr);
for i = 1:s
%    outarray(1+(i-1)*sr:sr*i)= sin((2*pi*scale(i))*t);%sin((2*pi*scale(i))*t);
     outarray(i,1:sr)= sin(2*pi*scale(i)*t)/sr/ns/s;%sin((2*pi*scale(i))*t);
end
figure()
imagesc(outarray)
title('Octavehigh')
%soundsc(outarray(8,1:2:sr/2),sr)

outarraymid = outarray;
%soundsc(outarraymid(8,1:2:sr/2),sr/2)
figure()
imagesc(outarraymid(:,1:sr/2))
title('Octavemiddle')

startnoteoffset = 0;%middleC: +3
isminor = 0;
sr = 44100;%per second
speed = 60;
ns =2;
s = 8;%
ds = 1;
majors=[0 2 4 5 7 9 11 12];%startoffset+3 : middle major C
minors=[0 2 3 5 7 8 10 12];
if(isminor == 1)
    scale = majors;
else
    scale = minors;
end
%the statement above is equivalent to 
sr = sr/ds;
t = 1:sr;
t = t/sr;
%unit = sr/ns*speed/60*s/8;
%fundamental 220Hz, Harmonics: 440Hz,  660Hz...
%Octave Example on a guitar:
%from 0 (2: 110)
%1 octave(1/2: 55), 2 octaves(4:440)
%1. E (full length; 66 cm)
%2. E' (1/2 of length; one octave up; 164 Hz; 33 cm; touching the string in the middle suppress the fundamental or first harmonic so you hear this pitch)
%3. B' (1/3 of length; an octave and a fifth up; 246 Hz; 22 cm; touching the string here suppresses the second harmonic so you hear this pitch)
%4. E'' (1/4 of length; two octaves up; 328 Hz; 16.5 cm; touching the string here suppresses the third harmonic so you hear this pitch)
%5. G''' (1/5 of length; two octaves and a third up; 410 Hz; 13.2 cm; touching the string here suppresses the fourth harmonic so you hear this pitch; in fact a tempered G''' is 415.305 Hz, so the tuner will register this note as "sharp" if you touch the string right above the fret)

startnote = 220*(2^(startnoteoffset+3));
scale = startnote * scale;
%note = ['8', '7', '8', '6', '3', '7', '5', '3', '-'; '8', '8', '7', '6', '8', '5', '7', '5'];
%note =['6', '-', '7', '5', '4', '3', '2', '1', '7', '1', '7', '-', '6', '5', '-'];
PLOT = 1;
%TTLS = MakeScale(1,startnoteoffset, isminor, sr, ds, s, '*');
line = 1;
Note = zeros(8,8);
Note(1,:) = [5,5,8,8,8,9,8,7];
Note(2,:) = [6,6,6,6,6,6,6,6];
Note(3,:) = [6,6,10,10,10,11,10,9];
Note(4,:) = [7,7,5,5,5,5,5,5];

Outarray = zeros(size(Note, 1), sr*s);
%outarray = zeros(1,sr*s);
outarray = zeros(s,sr);
outarraymid = zeros(s,sr);
for i = 1:s
%    outarray(1+(i-1)*sr:sr*i)= sin((2*pi*scale(i))*t);%sin((2*pi*scale(i))*t);
     outarray(i,1:sr)= sin(2*pi*scale(i)*t)/sr/ns/s;%sin((2*pi*scale(i))*t);
end
figure()
imagesc(outarray(:,1:8:sr/2))
title('Octavehigh')
%soundsc(outarray(8,1:8:sr/2),sr/2)

outarraymid = outarray;
figure()
imagesc(outarraymid(:,1:2:sr/2))
title('Octavemiddle')
%soundsc(outarraymid(8,1:2:sr),sr)

Note2 = zeros(8,8);
Note2(1,:) = [5,5,11,11,11,12,11,10];
Note2(2,:) = [9,9,6,6,6,6,6,6];
Note2(3,:) = [5,5,6,6,10,10,7,7];
Note2(4,:) = [8,8,8,8,8,8,8,8];
    %TTLS = MakeScale_1(startnoteoffset, isminor, PLOT, speed, sr, ns, s, note(2*line-1:2*line,:), unit)
    %TTLS = MakeScale(startnoteoffset, isminor, sr, ds, s, note)
    %
%Plot one cube(facelet)
if PLOT == 1
    figure(),
    display('Initialte the cube now(3*3*3).')
%        fprintf('Initialte the cube now(3*3*3).')
    d = 3
    R = rubgen(d,0);
    %R = RUBGEN_1(d,0);
    
    rubplot(R);
    display('Focus on the facelet facing you.')   
end

moves1 = {'x11','x31','y21','z11','z31','y14','y11','y11'};
figure(),
R1 = rubgen(d,0);
display('Forming letter K')

for i = 1:size(Note2,2)
    if Note(1,i) < 9
%        soundsc(outarraymid(mod(Note2(1,i),8),1:sr/2),sr/2)
    else
%        soundsc(outarray(mod(Note2(1,i),8),1:sr/2),sr)
    end 
    hold on
    pulsesep(1:pi/100:10*pi,1000);
    hold on
    if i <= length(moves1)
        R1 = rubrot(R1,moves1{i})
        rubplot(R1);
    end
    if i < size(Note2,2)
        hold on
    else
        hold off
    end
end
% R1 = rubgen(d,0);
% for l = 1:length(moves1)
%     R1 = rubrot(R1,moves1{i})
%     rubplot(R1);
% end

moves2 = {'y11','z11'};
figure(),
R2 = rubgen(d,0);
display('Forming letter r')

for i = 1:size(Note,2)
    if Note(1,i) < 9
%        soundsc(outarraymid(mod(Note(2,i),8),1:sr/2),sr/2)
    else
%        soundsc(outarray(mod(Note(2,i),8),1:sr/2),sr)
    end  
    hold on
    pulsesep(1:pi/100:10*pi,1000);
    hold on
    if i <= length(moves2)
        R2 = rubrot(R2,moves2{i})
        rubplot(R2);
    end
    if i < size(Note,2)
        hold on
    else
        hold off
    end
end

moves3 = {'z21','y11','y31'};
figure(),
R3 = rubgen(d,0);
display('Forming letter H')

for i = 1:size(Note,2)
    if Note(1,i) < 9
%        soundsc(outarraymid(mod(Note(3,i),8),1:sr/2),sr/2)
    else
%        soundsc(outarray(mod(Note(3,i),8),1:sr/2),sr)
    end  
    hold on
    pulsesep(1:pi/100:10*pi,1000);
    hold on
    if i <= length(moves3)
        R3 = rubrot(R3,moves3{i})
        rubplot(R3);
    end
    if i < size(Note,2)
        hold on
    else
        hold off
    end
end

moves4 = {'z11','z31','y11','y31','z11','z31','x31','x11'};
figure(),
R4 = rubgen(d,0);
display('Forming letter Q')

for i = 1:size(Note,2)
    if Note(1,i) < 9
%        soundsc(outarraymid(mod(Note(4,i),8),1:sr/2),sr/2)
    else
%        soundsc(outarray(mod(Note(4,i),8),1:sr/2),sr)
    end  
    hold on
    pulsesep(1:pi/100:10*pi,1000);
    hold on
    if i <= length(moves4)
        R4 = rubrot(R4,moves4{i})
        rubplot(R4);
    end
    if i < size(Note,2)
        hold on
    else
        hold off
    end
end

if PLOT == 1
    figure(),
    display('Initialte another cube now(3*3*3).')
%        fprintf('Initialte the cube now(3*3*3).')
    d = 3
    r = rubgen(d,0);
    %R = RUBGEN_1(d,0);
    
    rubplot(r);
    display('Focus on the facelet facing you.')   
end

moves1 = {'x11','x31','y21','z11','z31','y14','y11','y11'};
figure(),
r1 = rubgen(d,0);
display('Forming letter K')

for i = 1:size(Note2,2)
    if Note2(1,i) < 9
%        soundsc(outarraymid(mod(Note2(1,i),8),1:sr/2),sr/2)
    else
%        soundsc(outarray(mod(Note2(1,i),8),1:sr/2),sr)
    end 
    hold on
    pulsesep(1:pi/100:10*pi,1000);
    hold on
    if i <= length(moves1)
        r1 = rubrot(r1,moves1{i})
        rubplot(r1);
    end
    if i < size(Note2,2)
        hold on
    else
        hold off
    end
end
% R1 = rubgen(d,0);
% for l = 1:length(moves1)
%     R1 = rubrot(R1,moves1{i})
%     rubplot(R1);
% end

moves2 = {'y11','z11'};
figure(),
r2 = rubgen(d,0);
display('Forming letter r')

for i = 1:size(Note,2)
    if Note(1,i) < 9
%        soundsc(outarraymid(mod(Note2(2,i),8),1:sr/2),sr)
    else
%        soundsc(outarray(mod(Note2(2,i),8),1:sr/2),sr)
    end  
    hold on
    pulsesep(1:pi/100:10*pi,1000);
    hold on
    if i <= length(moves2)
        r2 = rubrot(r2,moves2{i})
        rubplot(r2);
    end
    if i < size(Note2,2)
        hold on
    else
        hold off
    end
end

moves3 = {'z21','y11','y31'};
figure(),
r3 = rubgen(d,0);
display('Forming letter H')

for i = 1:size(Note2,2)
    if Note2(1,i) < 9
%        soundsc(outarraymid(mod(Note2(3,i),8),1:sr/2),sr/2)
    else
%        soundsc(outarray(mod(Note2(3,i),8),1:sr/2),sr)
    end  
    hold on
    pulsesep(1:pi/100:10*pi,1000);
    hold on
    if i <= length(moves3)
        r3 = rubrot(r3,moves3{i})
        rubplot(r3);
    end
    if i < size(Note2,2)
        hold on
    else
        hold off
    end
end

moves4 = {'z11','z31','y11','y31','z11','z31','x31','x11'};
figure(),
r4 = rubgen(d,0);
display('Forming letter Q')

for i = 1:size(Note2,2)
    if Note2(1,i) < 9
%        soundsc(outarraymid(mod(Note2(4,i),8),1:sr/2),sr/2)
    else
%        soundsc(outarray(mod(Note2(4,i),8),1:sr/2),sr)
    end  
    hold on
    pulsesep(1:pi/100:10*pi,1000);
    hold on
    if i <= length(moves4)
        r4 = rubrot(R4,moves4{i})
        rubplot(r4);
    end
    if i < size(Note2,2)
        hold on
    else
        hold off
    end
end

for i = 1:size(Note,2)/2
    temp = mod(Note(i,:),8);
    for j = 1:8
%        soundsc(outarray(temp(j)+1,1:sr/2),sr/2);
        pulsesep(1:pi/100:10*pi,1000);
    end
end

for i = 1:size(Note2,2)/2
    temp = mod(Note2(i,:),8)+1;
    for j = 1:8
%        soundsc(outarray(temp(j),1:sr/2),sr/2);
        pulsesep(1:pi/100:10*pi,1000);
    end
end
%    X = [0 3 3 0;0 0 0 0;0 3 3 0;3 3 3 3; 0 3 3 0; 0 3 3 0];
%    Y = [0 0 3 3;0 3 3 0;0 0 3 3;0 3 3 0; 0 0 0 0; 3 3 3 3];
%    Z = [0 0 0 0;0 0 3 3;3 3 3 3;0 0 3 3; 0 0 3 3; 0 0 3 3];

    %Color = [hsv(1:ceil((length(hsv))/6):ceil(length(hsv)/6)*1);
    %hsv(1:ceil((length(hsv)/6):ceil(length(hsv)/6)*2);
    %hsv(1:ceil((length(hsv))/6):ceil(length(hsv)/6)*3;
    %hsv(1:ceil((length(hsv))/6):ceil(length(hsv)/6)*4);
    %hsv(1:ceil((length(hsvn))/6):ceil(length(hsv)/6)*5);
    %winter(1:ceil((length(hsv))/6):ceil(length(hsv)/6)*6)] Color
    %=[[cool(1) 0.5];[hot(1) 0.5];[spring(1) 0.5];[summer(1)
    %0.5];[autumn(1) 0.5]; [winter(1) 0.5]]
%     Color = [[spring(1) spring(1)]; [summer(1) summer(1)]; [autumn(1) autumn(1)]; [winter(1) winter(1)]];
%     figure,
%     fill3(X,Y,Z,Color')
% 
% 
%     for i = 1:size(Note, 1)
%         d = 3
%         R = rubgen(d,0);
%         for j = 1:size(Note, 2)
%  %           display('origin:')
%             NNote = Note(i, j);
%             if NNote == '-'
%                 assert(j>1);
%                 jj = j;
%                 while isempty(str2num(NNote))
%                     jj = jj-1;
%                 end
%                 Outarray(i,(jj-1)*sr+1:jj*sr) = Outarray(i,max((jj-2),0)*sr+1:(jj-1)*sr)
%             elseif NNote == '*'
%                 Outarray(i,(j-1)*sr+1:j*sr) = zeros(1,sr)
%                 move = 'xyz'
%                 Move = move;
%                 Move(1) = move(ceil(mod(NNote,8)/3));
%                 switch num2str(mod(str2num(NNote),3))
%                     case 0
%                         Move(2) = '1'
%                     case 1
%                         Move(2) = '2'
%                     otherwise
%                         Move(2) = '9'
%                     end
%                     Move(3) = '4';
%                     R = rubrot(R, move);
%                     figure(),
%                     rubplot(R);
%             elseif (NNote > 0 && NNote <= 8)
%    %         MakeScale_1(startnoteoffset, isminor, PLOT, speed, sr, ns, 1,
%    %         Note, unit)
%    
%    %         check= input(['input ', num2str((j+1)/2),'th note.'])
% %            subplot(3,3,j) PLOT = 1
%                 if NNote == 8
%                     Outarray(i,(j-1)*sr+1:j*sr) = outarray(1+(8-1)*sr:sr*i); 
%                 else 
%                     Outarray(i,(j-1)*sr+1:j*sr) = outarray(1+(NNote-1)*sr:sr*i); 
%                 end
%             elseif NNote > 8
%     %     %                fill3([1 2 2 1],[0 0 1 1],[0 0 0 0],TTLS(4, 1))
%     %                     pcolor([1 2; 1 2], [1 1; 2 2],
%     %                     reshape(Color(:,2), [2, 2]))
%                 Outarray(i,(j-1)*sr+1:j*sr) = NNote/8*outarray(1+(mod(NNote,8)-1)*sr:sr*i);
%             end
%             soundsc(Outarray(i, 1+(j-1)*sr:sr*j), unit, 8);
%             move  ='xyz';
%             Move = move;
%             Move(1) = move(ceil(mod(NNote,8)/3));
%             switch mod(NNote,3)
%                 case 0
%                     Move(2) = '1'
%                 case 1
%                     Move(2) = '2'
%                 otherwise
%                     Move(2) = '9'
%             end
%             Move(3) = '4'; 
%             R = rubrot(R, Move);
%             figure(),
%             rubplot(R);
%         end
%         figure(),
%         plot(Outarray(i,1:sr))
%         soundsc(Outarray(i,1:sr), unit, 8)
%     end
%     display('Now you can play the cube  by yourself.')          
%     E = GetEdges(R)
%     C = GetCorners(R);
%     digrub();
% 
% 
