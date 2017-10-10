function dx1dx1 = fse_dx1dx1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1)
%FSE_DX1DX1
%    DX1DX1 = FSE_DX1DX1(F1,F2,A1,A2,A3,A4,A7,ALPHA,BETA,X1)

%    This function was generated by the Symbolic Math Toolbox version 5.6.
%    28-Dec-2012 14:33:18

t367 = F2.*a1;
t365 = F1-t367;
t366 = alpha.*t365;
t368 = beta-t366;
t369 = exp(t368);
t370 = t369+1.0;
t371 = 1.0./t370;
t372 = a3-x1;
t373 = alpha.*t372;
t374 = 1.0./a3;
t375 = a7+1.0;
t376 = a2-1.0;
t377 = t365.*t371.*t376;
t378 = -beta-t373;
t379 = exp(t378);
t380 = t379+1.0;
t381 = 1.0./t380;
t385 = t374.*x1;
t382 = -t385+1.0;
t383 = 1.0./t375;
t384 = t382.^t383;
t386 = a3.*t375.*t384;
t387 = t377+t386;
t394 = alpha.*t387;
t388 = exp(-t394);
t389 = t388+1.0;
t390 = 1.0./t389;
t391 = t365.*t371.*t374.*t376.*t383;
t392 = t384+t391;
t393 = t392.^t375;
t395 = a3.*t390.*t393;
t396 = a3+t395-x1;
t397 = t371.*t381.*t396;
t398 = t377+t397;
t399 = -beta+t366;
t400 = exp(t399);
t401 = t400+1.0;
t402 = 1.0./t401;
t403 = alpha.*t398;
t404 = exp(t403);
t405 = t404+1.0;
t406 = 1.0./t405;
t407 = t383-1.0;
t408 = t382.^t407;
t409 = t392.^a7;
t410 = t390.*t408.*t409;
t411 = 1.0./t389.^2;
t412 = a3.*alpha.*t388.*t393.*t408.*t411;
t413 = t410+t412+1.0;
t414 = t371.*t381.*t413;
t415 = 1.0./t380.^2;
t416 = alpha.*t371.*t379.*t396.*t415;
t417 = t414+t416;
t418 = -beta+t373;
t419 = exp(t418);
t420 = t419+1.0;
t421 = t371.*t398.*t406;
t422 = -t377+t421+x1;
t423 = 1.0./t420;
t424 = a3.*t371.*t423;
t425 = t371.*t381.*t390.*t422;
t426 = t424+t425+x1;
t427 = t365.*t374.*t402.*x1;
t431 = a4.*t371.*t374.*t426.*(1.0./2.0);
t428 = t427-t431+x1;
t432 = alpha.*t428;
t429 = -beta-t432;
t430 = exp(t429);
t433 = t430+1.0;
t434 = t365.*t374.*t402;
t435 = 1.0./t405.^2;
t436 = alpha.*t371.*t398.*t404.*t417.*t435;
t447 = t371.*t406.*t417;
t437 = t436-t447+1.0;
t438 = t371.*t381.*t390.*t437;
t439 = 1.0./t420.^2;
t440 = a3.*alpha.*t371.*t419.*t439;
t448 = alpha.*t371.*t379.*t390.*t415.*t422;
t449 = alpha.*t371.*t381.*t388.*t408.*t411.*t422;
t441 = t438+t440-t448-t449+1.0;
t450 = a4.*t371.*t374.*t441.*(1.0./2.0);
t442 = t434-t450+1.0;
t443 = a3-t427+t431-x1;
t444 = alpha.*t443;
t445 = -beta+t444;
t446 = exp(t445);
dx1dx1 = t438+t440-t448-t449+t402.*(t442./t433+a3.*alpha.*t442.*t446.*1.0./(t446+1.0).^2+alpha.*t428.*t430.*1.0./t433.^2.*t442);
