t367 = nansum([F2 * a1, 0]);
t365 = F1-t367;
t366 = nansum([alpha * t365, 0]);
t368 = beta-t366;
t369 = exp(t368);
t370 = t369+1;
t371 = 1/t370;
t372 = a3-x1;
t373 = nansum([alpha * t372, 0]);
t374 = 1/a3;
t375 = a7+1;
t376 = a2-1;
t377_1 = nansum([t365 * t371, 0]);
t377 = nansum([t377_1 * t376,0]);
t378 = -beta-t373;
t379 = exp(t378);
t380 = t379+1;
t381 = 1/t380;
t385 = nansum([t374 * x1, 0]);
t382 = -t385+1;
t383 = 1/t375;
t384 = t382^t383;
t386_1 = nansum([a3 * t375, 0]);
t386 = nansum([t386_1 * t384,0]);
t387 = t377+t386;
t394 = nansum([alpha * t387, 0]);
t388 = exp(-t394);
t389 = t388+1;
t390 = 1/t389;
t391_1 = nansum([t365 * t371, 0]);
t391_2 = nansum([t391_1 * t374,0]);
t391_3 = nansum([t391_2 * t376,0]);
t391 = nansum([t391_3 * t383,0]);
t392 = t384+t391;
t393 = t392^t375;
t395_1 = nansum([a3 * t390, 0]);
t395 = nansum([t395_1 * t393,0]);
t396 = a3+t395-x1;
t397_1 = nansum([t371 * t381, 0]);
t397 = nansum([t397_1 * t396,0]);
t398 = t377+t397;
t399 = -beta+t366;
t400 = exp(t399);
t401 = t400+1;
t402 = 1/t401;
t403 = nansum([alpha * t398, 0]);
t404 = exp(t403);
t405 = t404+1;
t406 = 1/t405;
t407 = t383-1;
t408 = t382^t407;
t409 = t392^a7;
t410_1 = nansum([t390 * t408, 0]);
t410 = nansum([t410_1 * t409,0]);
t411 = 1/t389^2;
t412_1 = nansum([a3 * alpha, 0]);
t412_2 = nansum([t412_1 * t388,0]);
t412_3 = nansum([t412_2 * t393,0]);
t412_4 = nansum([t412_3 * t408,0]);
t412 = nansum([t412_4 * t411,0]);
t413 = t410+t412+1;
t414_1 = nansum([t371 * t381, 0]);
t414 = nansum([t414_1 * t413,0]);
t415 = 1/t380^2;
t416_1 = nansum([alpha * t371, 0]);
t416_2 = nansum([t416_1 * t379,0]);
t416_3 = nansum([t416_2 * t396,0]);
t416 = nansum([t416_3 * t415,0]);
t417 = t414+t416;
t418 = -beta+t373;
t419 = exp(t418);
t420 = t419+1;
t421_1 = nansum([t371 * t398, 0]);
t421 = nansum([t421_1 * t406,0]);
t422 = -t377+t421+x1;
t423 = 1/t420;
t424_1 = nansum([a3 * t371, 0]);
t424 = nansum([t424_1 * t423,0]);
t425_1 = nansum([t371 * t381, 0]);
t425_2 = nansum([t425_1 * t390,0]);
t425 = nansum([t425_2 * t422,0]);
t426 = t424+t425+x1;
t427_1 = nansum([t365 * t374, 0]);
t427_2 = nansum([t427_1 * t402,0]);
t427 = nansum([t427_2 * x1,0]);
t431_1 = nansum([a4 * t371, 0]);
t431_2 = nansum([t431_1 * t374,0]);
t431_3 = nansum([t431_2 * t426,0]);
t431 = nansum([t431_3 * (1/2),0]);
t428 = t427-t431+x1;
t432 = nansum([alpha * t428, 0]);
t429 = -beta-t432;
t430 = exp(t429);
t433 = t430+1;
t434_1 = nansum([t365 * t374, 0]);
t434 = nansum([t434_1 * t402,0]);
t435 = 1/t405^2;
t436_1 = nansum([alpha * t371, 0]);
t436_2 = nansum([t436_1 * t398,0]);
t436_3 = nansum([t436_2 * t404,0]);
t436_4 = nansum([t436_3 * t417,0]);
t436 = nansum([t436_4 * t435,0]);
t447_1 = nansum([t371 * t406, 0]);
t447 = nansum([t447_1 * t417,0]);
t437 = t436-t447+1;
t438_1 = nansum([t371 * t381, 0]);
t438_2 = nansum([t438_1 * t390,0]);
t438 = nansum([t438_2 * t437,0]);
t439 = 1/t420^2;
t440_1 = nansum([a3 * alpha, 0]);
t440_2 = nansum([t440_1 * t371,0]);
t440_3 = nansum([t440_2 * t419,0]);
t440 = nansum([t440_3 * t439,0]);
t448_1 = nansum([alpha * t371, 0]);
t448_2 = nansum([t448_1 * t379,0]);
t448_3 = nansum([t448_2 * t390,0]);
t448_4 = nansum([t448_3 * t415,0]);
t448 = nansum([t448_4 * t422,0]);
t449_1 = nansum([alpha * t371, 0]);
t449_2 = nansum([t449_1 * t381,0]);
t449_3 = nansum([t449_2 * t388,0]);
t449_4 = nansum([t449_3 * t408,0]);
t449_5 = nansum([t449_4 * t411,0]);
t449 = nansum([t449_5 * t422,0]);
t441 = t438+t440-t448-t449+1;
t450_1 = nansum([a4 * t371, 0]);
t450_2 = nansum([t450_1 * t374,0]);
t450_3 = nansum([t450_2 * t441,0]);
t450 = nansum([t450_3 * (1/2),0]);
t442 = t434-t450+1;
t443 = a3-t427+t431-x1;
t444 = nansum([alpha * t443, 0]);
t445 = -beta+t444;
t446 = exp(t445);
M2_1 = nansum([alpha * t428, 0]);
M2_2 = nansum([M2_1 * t430,0]);
M2_3 = nansum([M2_2 * (1/t433^2),0]);
M2 = nansum([M2_3 * t442,0]);
M1_1 = nansum([a3 * alpha, 0]);
M1_2 = nansum([M1_1 * t442,0]);
M1_3 = nansum([M1_2 * t446,0]);
M1 = nansum([M1_3 * (1/(t446+1)^2),0]);
M = (t442/t433)+M1+M2;
MX = nansum([t402 * M, 0]);
dx1dx1 = t438+t440-t448-t449+MX;