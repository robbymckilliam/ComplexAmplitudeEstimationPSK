prologues := 2;
filenametemplate "%j-%c.mps";

input TEX; 
TEXPRE("%&latex" & char(10) & 
"\documentclass{minimal}\begin{document}"); 
TEXPOST("\end{document}"); 

input graph;
input sarith;
picture ORIGINAL_PICTURE;
ORIGINAL_PICTURE := currentpicture;

picture dd; dd :=
dashpattern(on 4bp off 2bp on 1bp off 2bp);

picture diamond;
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm rotated 45;
diamond := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture square;
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm;
square := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture circle; 
draw fullcircle scaled 0.19cm;
circle := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture cross; 
draw (-0.09cm,0.0)--(0.09cm,0.0);
draw (0,-0.09cm)--(0,0.09cm);
cross := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture dot; 
fill fullcircle scaled 0.15cm;
dot := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture times; times := cross rotated 45;


%%%%%%%%%%%%%
%%%%%Legend picture

picture legend;
sx = 0.2cm;
sh = 0.2cm;
ups = 0.35cm;

ls:=0;
label.rt(btex $Mackenthun$ etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw dot shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=1;
label.rt(btex $Unmodulated Carrier$ etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw cross shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=2;
label.rt(btex $Viterbi Viterbi$ etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

legend := currentpicture; currentpicture := ORIGINAL_PICTURE;


plotheight = 16cm;
plotwidth = 11cm;

legendx = 0.1cm;
legendy = 0.1cm;

M = 8;
ymax = (3.3/M)*(3.3/M);

%%%%%%%%%%%%
%%%% Begin plot of amplitude error

beginfig(1);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (30cm,1.3cm);
  glabel.lft(btex MSE amplitude etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  %setrange(whatever,"1.1e-6", whatever,90);
  %setrange(4.9,"1e-13", 22,2);

%  gdraw (0, "1e-13")--(0, 2) dashed evenly; %glabel.bot(btex AWGN channel capacity etex, 0) rotated 90;

  gdraw "MackenthunNonCoherentM8L1000a"; gdraw "MackenthunNonCoherentM8L1000a" plot(dot); glabel.llft(btex $L = 1000$ etex, 1000);
  gdraw "MackenthunNonCoherentM8L100a"; gdraw "MackenthunNonCoherentM8L100a" plot(dot); glabel.llft(btex $L = 100$ etex, 1000);
  gdraw "MackenthunNonCoherentM8L10a"; gdraw "MackenthunNonCoherentM8L10a" plot(dot); glabel.llft(btex $L = 10$ etex, 1000);
  %gdraw "MackenthunNonCoherentM8L10000absP1000a"; gdraw "MackenthunNonCoherentM8L10000absP1000a" plot(dot); glabel.llft(btex $L = 10000$ etex, 1000);

  gdraw "ViterbiViterbiM8L1000a"; gdraw "ViterbiViterbiM8L1000a" plot(circle); 
  gdraw "ViterbiViterbiM8L100a"; gdraw "ViterbiViterbiM8L100a" plot(circle); 
  gdraw "ViterbiViterbiM8L10a"; gdraw "ViterbiViterbiM8L10a" plot(circle);

  gdraw "MackenthunCoherentM8L1000absP1000a"; 
  gdraw "MackenthunCoherentM8L100absP100a"; 
  gdraw "MackenthunCoherentM8L10absP10a"; 


draw legend shifted(legendx,legendy);

endgraph;

clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;


%%%%%%%%%%%%
%%%% Begin plot of phase error

beginfig(2);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0cm,1.3cm);
  glabel.lft(btex MSE phase etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  setrange(whatever,whatever, whatever,ymax);
  %setrange(4.9,"1e-13", 22,2);

%  gdraw (0, "1e-13")--(0, 2) dashed evenly; %glabel.bot(btex AWGN channel capacity etex, 0) rotated 90;

  gdraw "MackenthunNonCoherentM8L1000p"; gdraw "MackenthunNonCoherentM8L1000p" plot(dot); glabel.llft(btex $L = 1000$ etex, 1000);
  gdraw "MackenthunNonCoherentM8L100p"; gdraw "MackenthunNonCoherentM8L100p" plot(dot); glabel.llft(btex $L = 100$ etex, 1000);
  gdraw "MackenthunNonCoherentM8L10p"; gdraw "MackenthunNonCoherentM8L10p" plot(dot); glabel.llft(btex $L = 10$ etex, 1000);
  %gdraw "MackenthunNonCoherentM8L10000absP1000p"; gdraw "MackenthunNonCoherentM8L10000absP1000p" plot(dot); glabel.llft(btex $L = 10000$ etex, 1000);

   gdraw "ViterbiViterbiM8L1000p"; gdraw "ViterbiViterbiM8L1000p" plot(circle); 
  gdraw "ViterbiViterbiM8L100p"; gdraw "ViterbiViterbiM8L100p" plot(circle); 
  gdraw "ViterbiViterbiM8L10p"; gdraw "ViterbiViterbiM8L10p" plot(circle);

  gdraw "MackenthunCoherentM8L1000absP1000p"; 
  gdraw "MackenthunCoherentM8L100absP100p"; 
  gdraw "MackenthunCoherentM8L10absP10p";

draw legend shifted(legendx,legendy);

endgraph;

clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;


%%%%%%%%%%%%
%%%% Begin plot of complex error

beginfig(3);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0cm,1.3cm);
  glabel.lft(btex MSE complex gain etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  %setrange(whatever,"1.1e-6", whatever,90);
  %setrange(4.9,"1e-13", 22,2);

%  gdraw (0, "1e-13")--(0, 2) dashed evenly; %glabel.bot(btex AWGN channel capacity etex, 0) rotated 90;

  gdraw "MackenthunNonCoherentM8L1000c"; gdraw "MackenthunNonCoherentM8L1000c" plot(dot); glabel.llft(btex $L = 1000$ etex, 1000);
  gdraw "MackenthunNonCoherentM8L100c"; gdraw "MackenthunNonCoherentM8L100c" plot(dot); glabel.llft(btex $L = 100$ etex, 1000);
  gdraw "MackenthunNonCoherentM8L10c"; gdraw "MackenthunNonCoherentM8L10c" plot(dot); glabel.llft(btex $L = 10$ etex, 1000);
  %gdraw "MackenthunNonCoherentM8L10000absP1000c"; gdraw "MackenthunNonCoherentM8L10000absP1000c" plot(dot); glabel.llft(btex $L = 10000$ etex, 1000);

  gdraw "ViterbiViterbiM8L1000c"; gdraw "ViterbiViterbiM8L1000c" plot(circle); 
  gdraw "ViterbiViterbiM8L100c"; gdraw "ViterbiViterbiM8L100c" plot(circle); 
  gdraw "ViterbiViterbiM8L10c"; gdraw "ViterbiViterbiM8L10c" plot(circle);
  
  gdraw "MackenthunCoherentM8L1000absP1000c"; 
  gdraw "MackenthunCoherentM8L100absP100c"; 
  gdraw "MackenthunCoherentM8L10absP10c";

draw legend shifted(legendx,legendy);

endgraph;

clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;

end;