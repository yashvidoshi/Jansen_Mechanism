clear; close all; clc;

if exist("fsolve", "file") == 2
    fprintf("Solver: MATLAB fsolve from Optimization Toolbox.\n");
else
    fprintf("Solver: included damped Newton solver because fsolve is unavailable.\n");
end

%% 1. Parameters
% Appendix nominal set: [11,45,36,34,48.5,41.5,60.5,41.5,42,43,26.5,54.5].
L = struct();
L.names = ["L1","L2","L3","L4","L5","L6","L7","L8","L9","L10","L11","L12"];
L.val = [11.29, 45.0, 36.0, 32.93, 48.5, 41.5, 60.5, 41.78, 42.0, 43.0, 26.5, 54.5].';

sim = struct();
sim.N = 361;                    % one full crank turn
sim.omega = 2*pi/2.0;           % rad/s; 2 s per gait cycle
sim.phaseOffset = 0.0;          % radians; changes starting point only
sim.gaitFrameRotationDeg = -11.4; % rotates mechanism frame to gait/world frame
sim.showPaperFigure4 = true;
sim.showNinePatternValidation = true;
sim.runAnimation = true;
sim.animationFrameStep = 3;      % higher value makes animation faster
sim.showDiagnosticPlots = false; % set true for separate trajectory/velocity/residual figures
sim.showMechanismFrames = false;
sim.NAdvanced = 91;             % lower resolution for sensitivity sweeps
sim.runSensitivity = false;      % set true for adjustable-link sensitivity plots
sim.runLeastSquaresDemo = false; % set true for span-to-link least-squares demo

fprintf("Modified Jansen gait trainer simulation\n");
fprintf("Using adjustable values: L1 = %.2f cm, L4 = %.2f cm, L8 = %.2f cm\n", ...
    L.val(1), L.val(4), L.val(8));

%% 2. Main one-cycle simulation
main = simulateJansenCycle(L.val, sim.N, sim.omega, sim.phaseOffset, [], ...
    sim.gaitFrameRotationDeg*pi/180);

spanX = spanOf(main.PE(1,:));
spanY = spanOf(main.PE(2,:));
areaPE = abs(polyarea(main.PE(1,:), main.PE(2,:)));
fprintf("End-effector x-span = %.2f cm\n", spanX);
fprintf("End-effector y-span = %.2f cm\n", spanY);
fprintf("Closed-loop area     = %.2f cm^2\n", areaPE);
fprintf("Max loop residual    = %.3e cm\n\n", max(main.resnorm));

%% 3. Paper-style reference curve for visual comparison
% The actual human marker data used in the papers is not inside the PDFs.
% This pchip curve recreates the reported gait envelope:
% desired span [xspan, yspan] = [50.02, 12.81] cm.
targetSpan = [50.02; 12.81];
ref = makeReferenceGait(main.gaitCycle, targetSpan);
refAligned = alignCurveToSimulation(ref, main.PE);

%% 4. Paper-style outputs
if sim.showPaperFigure4
    plotPaperStyleFigure4(main, refAligned, L);
end

if sim.showNinePatternValidation
    validation = plotNinePatternValidation(L.val, sim);
end

if sim.runAnimation
    animateJansenMechanism(main, L, sim.animationFrameStep);
end

%% 5. Optional diagnostic plots
if sim.showDiagnosticPlots
    plotEndEffector(main, refAligned, targetSpan, L);
    plotGaitCurves(main, refAligned);
    plotVelocityCurves(main);
    plotConstraintResiduals(main);
end

if sim.showMechanismFrames
    plotMechanismSnapshots(main, L);
end

%% 6. Optional advanced plots: sensitivity and least-squares link-span mapping
if sim.runSensitivity
    sensitivity = runAdjustableLinkSensitivity(L.val, sim);
    plotSensitivity(sensitivity);
end

if sim.runLeastSquaresDemo
    lsDemo = runLeastSquaresSpanMap(L.val, sim, targetSpan);
    plotLeastSquaresDemo(lsDemo, targetSpan);
end

%% 7. Data exported to the MATLAB workspace
result = struct();
result.lengths_cm = L;
result.main = main;
result.reference = refAligned;
result.targetSpan_cm = targetSpan;
if exist("validation", "var")
    result.validation = validation;
end
if exist("sensitivity", "var")
    result.sensitivity = sensitivity;
end
if exist("lsDemo", "var")
    result.leastSquaresDemo = lsDemo;
end

disp("Done. The struct named 'result' contains trajectories, angles, spans, and residuals.");

%% Local functions

function out = simulateJansenCycle(L, N, omega, phaseOffset, qStart, frameRotation)
    if nargin < 6
        frameRotation = 0;
    end

    crank = linspace(0, 2*pi, N);
    t = crank / omega;
    gaitCycle = 100 * crank / (2*pi);

    solver = makeLoopSolver();

    theta = nan(12, N);
    qHist = nan(10, N);
    exitflag = nan(1, N);
    resnorm = nan(1, N);

    P0 = nan(2, N); P1 = P0; P2 = P0; P3 = P0; P4 = P0;
    P5 = P0; P6 = P0; PE = P0;

    q = qStart;
    for k = 1:N
        theta1 = crank(k) + phaseOffset;

        if isempty(q) || any(~isfinite(q))
            q = geometricInitialGuess(L, theta1);
        end

        loopFun = @(qq) loopEquations(qq, theta1, L);
        loopJac = @(qq) loopJacobian(qq, L);
        [qCandidate, fval, flag] = solveLoopAngles(loopFun, loopJac, q, solver);
        if flag <= 0 || norm(fval) > 1e-6
            qGeom = geometricInitialGuess(L, theta1);
            [qCandidate2, fval2, flag2] = solveLoopAngles(loopFun, loopJac, qGeom, solver);
            if norm(fval2) < norm(fval)
                qCandidate = qCandidate2;
                fval = fval2;
                flag = flag2;
            end
        end

        q = unwrapNear(qCandidate(:), q(:));
        qHist(:,k) = q;
        exitflag(k) = flag;
        resnorm(k) = norm(fval);

        th = nan(12,1);
        th(1) = theta1;
        th(2) = q(1);
        th(3) = q(2);
        th(4) = 0;
        th(5) = q(3);
        th(6) = q(4);
        th(7) = q(5);
        th(8) = q(6);
        th(9) = q(7);
        th(10) = q(8);
        th(11) = q(9);
        th(12) = q(10);
        theta(:,k) = th;

        pos = positionsFromAngles(L, th);
        P0(:,k) = pos.P0; P1(:,k) = pos.P1; P2(:,k) = pos.P2;
        P3(:,k) = pos.P3; P4(:,k) = pos.P4; P5(:,k) = pos.P5;
        P6(:,k) = pos.P6; PE(:,k) = pos.PE;
    end

    raw = struct("P0", P0, "P1", P1, "P2", P2, "P3", P3, "P4", P4, ...
        "P5", P5, "P6", P6, "PE", PE);

    if abs(frameRotation) > 0
        R = [cos(frameRotation), -sin(frameRotation); ...
             sin(frameRotation),  cos(frameRotation)];
        P0 = R*P0; P1 = R*P1; P2 = R*P2; P3 = R*P3; P4 = R*P4;
        P5 = R*P5; P6 = R*P6; PE = R*PE;
    end

    vx = gradient(PE(1,:), t);
    vy = gradient(PE(2,:), t);
    speed = hypot(vx, vy);

    out = struct();
    out.t = t;
    out.crank = crank;
    out.gaitCycle = gaitCycle;
    out.theta = theta;
    out.qHist = qHist;
    out.exitflag = exitflag;
    out.resnorm = resnorm;
    out.P0 = P0; out.P1 = P1; out.P2 = P2; out.P3 = P3; out.P4 = P4;
    out.P5 = P5; out.P6 = P6; out.PE = PE;
    out.rawMechanismFrame = raw;
    out.frameRotation_rad = frameRotation;
    out.vx = vx; out.vy = vy; out.speed = speed;
    out.span = [spanOf(PE(1,:)); spanOf(PE(2,:))];
    out.area = abs(polyarea(PE(1,:), PE(2,:)));
end

function F = loopEquations(q, theta1, L)
%LOOPEQUATIONS Five vector-loop closures, two scalar equations each.
    u = @(a) [cos(a); sin(a)];

    theta2 = q(1); theta3 = q(2); theta5 = q(3); theta6 = q(4);
    theta7 = q(5); theta8 = q(6); theta9 = q(7); theta10 = q(8);
    theta11 = q(9); theta12 = q(10);

    e1 = u(theta1);  e2 = u(theta2);  e3 = u(theta3);  e4 = [1;0];
    e5 = u(theta5);  e6 = u(theta6);  e7 = u(theta7);  e8 = u(theta8);
    e9 = u(theta9);  e10 = u(theta10); e11 = u(theta11); e12 = u(theta12);

    Fupper = L(1)*e1 + L(2)*e2 - L(3)*e3 - L(4)*e4;
    Flower = L(1)*e1 + L(7)*e7 - L(8)*e8 - L(4)*e4;
    Ftriangle = L(3)*e3 + L(5)*e5 - L(6)*e6;
    Fpara = L(8)*e8 + L(9)*e9 - L(6)*e6 - L(10)*e10;
    Ffoot = L(11)*e11 + L(12)*e12 - L(9)*e9;

    F = [Fupper; Flower; Ftriangle; Fpara; Ffoot];
end

function J = loopJacobian(q, L)
%LOOPJACOBIAN Analytic Jacobian dF/dq for the five vector-loop equations.
    d = @(a) [-sin(a); cos(a)];

    theta2 = q(1); theta3 = q(2); theta5 = q(3); theta6 = q(4);
    theta7 = q(5); theta8 = q(6); theta9 = q(7); theta10 = q(8);
    theta11 = q(9); theta12 = q(10);

    J = zeros(10,10);

    % Upper loop: L1*e1 + L2*e2 - L3*e3 - L4*e4 = 0.
    J(1:2,1) =  L(2)*d(theta2);
    J(1:2,2) = -L(3)*d(theta3);

    % Lower loop: L1*e1 + L7*e7 - L8*e8 - L4*e4 = 0.
    J(3:4,5) =  L(7)*d(theta7);
    J(3:4,6) = -L(8)*d(theta8);

    % Coupler triangle: L3*e3 + L5*e5 - L6*e6 = 0.
    J(5:6,2) =  L(3)*d(theta3);
    J(5:6,3) =  L(5)*d(theta5);
    J(5:6,4) = -L(6)*d(theta6);

    % Parallelogram-like loop: L8*e8 + L9*e9 - L6*e6 - L10*e10 = 0.
    J(7:8,4) = -L(6)*d(theta6);
    J(7:8,6) =  L(8)*d(theta8);
    J(7:8,7) =  L(9)*d(theta9);
    J(7:8,8) = -L(10)*d(theta10);

    % Foot triangle: L11*e11 + L12*e12 - L9*e9 = 0.
    J(9:10,7) = -L(9)*d(theta9);
    J(9:10,9) =  L(11)*d(theta11);
    J(9:10,10) = L(12)*d(theta12);
end

function solver = makeLoopSolver()
%MAKELOOPSOLVER Use fsolve when available; otherwise use a built-in solver.
    solver = struct();
    solver.useFsolve = exist("fsolve", "file") == 2;
    solver.tolF = 1e-10;
    solver.tolStep = 1e-11;
    solver.maxIter = 80;
    solver.maxLineSearch = 16;

    if solver.useFsolve
        solver.options = optimoptions("fsolve", ...
            "Display", "off", ...
            "FunctionTolerance", solver.tolF, ...
            "StepTolerance", solver.tolStep, ...
            "OptimalityTolerance", solver.tolF, ...
            "MaxIterations", 120, ...
            "MaxFunctionEvaluations", 1500);
    else
        solver.options = [];
    end
end

function [x, F, exitflag] = solveLoopAngles(fun, jac, x0, solver)
%SOLVELOOPANGLES Nonlinear solve wrapper.
    if solver.useFsolve
        [x, F, exitflag] = fsolve(fun, x0, solver.options);
        return;
    end

    [x, F, exitflag] = dampedNewtonSolve(fun, jac, x0, solver);
end

function [x, F, exitflag] = dampedNewtonSolve(fun, jac, x0, solver)
%DAMPEDNEWTONSOLVE Small square-system Newton solver for this linkage.
% It is included so the script remains runnable without Optimization Toolbox.
    x = x0(:);
    F = fun(x);
    nrm = norm(F);
    exitflag = 0;

    for iter = 1:solver.maxIter
        if nrm < solver.tolF
            exitflag = 1;
            return;
        end

        J = jac(x);
        step = -J \ F;
        if any(~isfinite(step)) || norm(step) > 5
            step = -pinv(J) * F;
        end

        alpha = 1.0;
        accepted = false;
        for ls = 1:solver.maxLineSearch
            xt = x + alpha*step;
            Ft = fun(xt);
            nt = norm(Ft);
            if nt < nrm || nt < solver.tolF
                x = xt;
                F = Ft;
                nrm = nt;
                accepted = true;
                break;
            end
            alpha = 0.5*alpha;
        end

        if ~accepted
            x = x + alpha*step;
            F = fun(x);
            nrm = norm(F);
        end

        if norm(alpha*step) < solver.tolStep*(1 + norm(x))
            exitflag = double(nrm < 1e-7);
            return;
        end
    end
end

function pos = positionsFromAngles(L, th)
%POSITIONSFROMANGLES Forward kinematics from solved link orientations.
    u = @(a) [cos(a); sin(a)];

    P0 = [0;0];
    P3 = [L(4);0];
    P1 = P0 + L(1)*u(th(1));
    P2 = P1 + L(2)*u(th(2));
    P5 = P1 + L(7)*u(th(7));
    P4 = P2 + L(5)*u(th(5));
    P6 = P5 + L(9)*u(th(9));
    PE = P5 + L(11)*u(th(11));

    pos = struct("P0", P0, "P1", P1, "P2", P2, "P3", P3, ...
        "P4", P4, "P5", P5, "P6", P6, "PE", PE);
end

function q = geometricInitialGuess(L, theta1)
%GEOMETRICINITIALGUESS Build the physical assembly mode from circle intersections.
% This is not the solver; it selects the same branch before fsolve refines
% the vector-loop equations.
    P0 = [0;0];
    P3 = [L(4);0];
    P1 = P0 + L(1)*[cos(theta1); sin(theta1)];

    P2c = circleIntersections(P1, L(2), P3, L(3));
    P2 = pickBy(P2c, "maxY");

    P5c = circleIntersections(P1, L(7), P3, L(8));
    P5 = pickBy(P5c, "minY");

    P4c = circleIntersections(P2, L(5), P3, L(6));
    P4 = pickBy(P4c, "maxX");

    P6c = circleIntersections(P5, L(9), P4, L(10));
    P6 = pickBy(P6c, "maxX");
    if P6(2) > min(P4(2), P5(2))
        P6 = pickBy(P6c, "minY");
    end

    PEc = circleIntersections(P5, L(11), P6, L(12));
    PE = pickBy(PEc, "minY");

    ang = @(a,b) atan2(b(2)-a(2), b(1)-a(1));
    q = [
        ang(P1, P2)
        ang(P3, P2)
        ang(P2, P4)
        ang(P3, P4)
        ang(P1, P5)
        ang(P3, P5)
        ang(P5, P6)
        ang(P4, P6)
        ang(P5, PE)
        ang(PE, P6)
    ];
end

function pts = circleIntersections(c1, r1, c2, r2)
%CIRCLEINTERSECTIONS Return the two intersection points of two circles.
    dvec = c2 - c1;
    d = norm(dvec);
    if d < eps
        error("Coincident circle centers in initial guess.");
    end
    if d > r1 + r2 || d < abs(r1 - r2)
        error("The selected link lengths cannot assemble at this crank angle.");
    end

    a = (r1^2 - r2^2 + d^2) / (2*d);
    h2 = max(r1^2 - a^2, 0);
    h = sqrt(h2);
    ex = dvec / d;
    ey = [-ex(2); ex(1)];
    p = c1 + a*ex;
    pts = [p + h*ey, p - h*ey];
end

function p = pickBy(pts, mode)
%PICKBY Select one of two branch points by geometric criterion.
    switch mode
        case "maxY"
            [~, idx] = max(pts(2,:));
        case "minY"
            [~, idx] = min(pts(2,:));
        case "maxX"
            [~, idx] = max(pts(1,:));
        case "minX"
            [~, idx] = min(pts(1,:));
        otherwise
            error("Unknown branch-selection mode.");
    end
    p = pts(:,idx);
end

function q = unwrapNear(qNew, qOld)
%UNWRAPNEAR Keep periodic angle solutions close to the previous step.
    if isempty(qOld) || any(~isfinite(qOld))
        q = qNew;
        return;
    end
    q = qNew + 2*pi*round((qOld - qNew)/(2*pi));
end

function ref = makeReferenceGait(gaitCycle, span)
%MAKEREFERENCEGAIT Smooth gait-like reference with the paper-reported spans.
% This is for visual validation only; it is not the unpublished human dataset.
    s = gaitCycle(:).' / 100;
    knot = [0.00 0.08 0.16 0.28 0.42 0.58 0.72 0.86 1.00];
    xShape = [0.48 0.42 0.25 -0.15 -0.50 -0.35 -0.10 0.20 0.48];
    yShape = [0.55 0.95 1.00 0.35 0.15 0.06 0.00 0.18 0.55];
    x = pchip(knot, xShape, s);
    y = pchip(knot, yShape, s);

    x = span(1) * (x - min(x)) / spanOf(x);
    y = span(2) * (y - min(y)) / spanOf(y);
    x = x - mean(x);
    y = y - min(y);

    ref = [x; y];
end

function refAligned = alignCurveToSimulation(ref, PE)
%ALIGNCURVETOSIMULATION Translate reference to same lower-left envelope.
    refAligned = ref;
    refAligned(1,:) = ref(1,:) - mean(ref(1,:)) + mean(PE(1,:));
    refAligned(2,:) = ref(2,:) - min(ref(2,:)) + min(PE(2,:));
end

function plotPaperStyleFigure4(main, ref, L)
%PLOTPAPERSTYLEFIGURE4 Match the composite mechanism/gait plot in the paper.
    fig = figure("Name", "Paper-style Fig. 4: mechanism and gait curves", ...
        "Color", "w", "Position", [80 80 1120 560]);

    axMech = axes(fig, "Position", [0.06 0.16 0.53 0.76]);
    axes(axMech);
    k = 22;
    drawPaperMechanism(main, k, L);
    hold on;

    markerIdx = 1:5:numel(main.gaitCycle);
    plot(ref(1,markerIdx), ref(2,markerIdx), "k*", ...
        "MarkerSize", 5.5, "LineWidth", 1.1);
    plot(main.PE(1,markerIdx), main.PE(2,markerIdx), "g", ...
        "MarkerSize", 4.0, "LineWidth", 1.1);
    text(main.PE(1,end)+1.0, main.PE(2,end), "Endpoint", ...
        "FontSize", 10, "VerticalAlignment", "middle");

    lengthText = makeLengthListText(L.val);
    xText = min([main.P0(1,:), main.PE(1,:)]) - 6;
    yText = max([main.P2(2,:), main.P4(2,:)]) - 4;
    text(xText, yText, lengthText, "FontSize", 9, "FontName", "Consolas", ...
        "VerticalAlignment", "top", "Interpreter", "none");

    axis equal;
    axis off;
    title("Parameterized 12-link modified Jansen mechanism", ...
        "FontSize", 11, "FontWeight", "normal");

    axX = axes(fig, "Position", [0.67 0.58 0.28 0.30]);
    predX = main.PE(1,:) - mean(main.PE(1,:));
    refX = ref(1,:) - mean(ref(1,:));
    plot(main.gaitCycle, refX, "m--", "LineWidth", 1.5); hold on;
    plot(main.gaitCycle, predX, "k-", "LineWidth", 1.8);
    grid on;
    xlim([0 100]);
    ylim([-40 40]);
    ylabel("x-axis (cm)");
    set(gca, "FontSize", 9);

    axY = axes(fig, "Position", [0.67 0.20 0.28 0.30]);
    predY = main.PE(2,:) - min(main.PE(2,:));
    refY = ref(2,:) - min(ref(2,:));
    plot(main.gaitCycle, refY, "b--", "LineWidth", 1.5); hold on;
    plot(main.gaitCycle, predY, "r-", "LineWidth", 1.8);
    grid on;
    xlim([0 100]);
    ylim([0 15]);
    xlabel("Gait Cycle (%)");
    ylabel("y-axis (cm)");
    legend("meta-trajectory", "predicted trajectory", ...
        "Location", "northeast", "FontSize", 8);
    set(gca, "FontSize", 9);

    caption = sprintf(['Fig. 4  Optimized to the gait envelope, the structure produces an endpoint trajectory with ' ...
        'x-span %.2f cm and y-span %.2f cm. Crosses denote the reference envelope; circles/solid curves denote the ' ...
        'simulated mechanism trajectory under constant crank speed.'], main.span(1), main.span(2));
    annotation(fig, "textbox", [0.08 0.015 0.86 0.10], ...
        "String", caption, ...
        "EdgeColor", "none", "FontSize", 10, "FontWeight", "bold");

    axes(axMech);
end

function drawPaperMechanism(main, k, L)
%DRAWPAPERMECHANISM Draw one pose with labels similar to the paper figure.
    linePairs = {
        "P0","P1","L_1"; "P1","P2","L_2"; "P3","P2","L_3"; "P0","P3","L_4";
        "P2","P4","L_5"; "P3","P4","L_6"; "P1","P5","L_7"; "P3","P5","L_8";
        "P5","P6","L_9"; "P4","P6","L_{10}"; "P5","PE","L_{11}"; "PE","P6","L_{12}"};

    for i = 1:size(linePairs,1)
        A = main.(linePairs{i,1})(:,k);
        B = main.(linePairs{i,2})(:,k);
        plot([A(1), B(1)], [A(2), B(2)], "-", "LineWidth", 1.5);
        hold on;
        mid = 0.52*A + 0.48*B;
        text(mid(1), mid(2), linePairs{i,3}, ...
            "FontSize", 11, "FontWeight", "bold", "Interpreter", "tex");
    end

    pointNames = ["P0","P1","P2","P3","P4","P5","P6","PE"];
    P = zeros(2, numel(pointNames));
    for i = 1:numel(pointNames)
        P(:,i) = main.(pointNames(i))(:,k);
    end
    plot(P(1,:), P(2,:), "ks", "MarkerFaceColor", "y", "MarkerSize", 6);

    xMargin = 8;
    yMargin = 8;
    xlim([min([P(1,:), main.PE(1,:)])-xMargin, max([P(1,:), main.PE(1,:)])+xMargin]);
    ylim([min([P(2,:), main.PE(2,:)])-yMargin, max([P(2,:), main.PE(2,:)])+yMargin]);

    % Put link labels a little away from the length-list text.
    unused = L; %#ok<NASGU>
end

function txt = makeLengthListText(L)
%MAKELENGTHLISTTEXT Link length annotation used in the paper-style figure.
    lines = strings(14,1);
    lines(1) = "Unit: cm";
    for i = 1:12
        lines(i+1) = sprintf("L%-2d = %4.1f", i, L(i));
    end
    txt = strjoin(lines, newline);
end

function validation = plotNinePatternValidation(Lnom, sim)
%PLOTNINEPATTERNVALIDATION Create a 3x3 paper-style RMSE validation grid.
% The PDFs do not contain the original 113-subject ankle database. This grid
% uses nine gait-envelope variants and solves the actual mechanism for each.
    xSpanGrid = [56 53 50; 53 50 47; 50 47 44];
    ySpanGrid = [14.2 13.6 13.0; 13.6 13.0 12.4; 13.0 12.4 11.8];
    L1Grid = Lnom(1) + [0.45 0.25 0.05; 0.25 0.00 -0.20; 0.05 -0.20 -0.45];
    L4Grid = Lnom(4) + [-0.90 -0.50 -0.10; -0.50 0.00 0.50; -0.10 0.50 0.90];
    L8Grid = Lnom(8) + [0.90 0.50 0.10; 0.50 0.00 -0.50; 0.10 -0.50 -0.90];

    fig = figure("Name", "Paper-style 3x3 validation grid", ...
        "Color", "w", "Position", [120 90 1040 620]);
    tl = tiledlayout(3,3, "Padding", "compact", "TileSpacing", "compact");

    validation = struct();
    validation.patterns = repmat(struct("L", [], "reference", [], ...
        "simulation", [], "rmse", [], "span", []), 3, 3);

    for row = 1:3
        for col = 1:3
            L = Lnom;
            L(1) = L1Grid(row,col);
            L(4) = L4Grid(row,col);
            L(8) = L8Grid(row,col);
            out = simulateJansenCycle(L, sim.NAdvanced, sim.omega, sim.phaseOffset, [], ...
                sim.gaitFrameRotationDeg*pi/180);

            ref = makeReferenceGait(out.gaitCycle, [xSpanGrid(row,col); ySpanGrid(row,col)]);
            ref = alignCurveToSimulation(ref, out.PE);

            [simPanel, refPanel] = panelNormalizeCurves(out.PE, ref);
            rmse = sqrt(mean(sum((simPanel - refPanel).^2, 1)));

            ax = nexttile;
            markerIdx = 1:2:numel(out.gaitCycle);
            plot(refPanel(1,markerIdx), refPanel(2,markerIdx), "k+", ...
                "MarkerSize", 5.0, "LineWidth", 1.0); hold on;
            plot(simPanel(1,markerIdx), simPanel(2,markerIdx), "gs", ...
                "MarkerSize", 4.0, "LineWidth", 1.0);
            grid on;
            xlim([0 70]);
            ylim([0 20]);
            text(5, 17, sprintf("RMSE=%.2f", rmse), ...
                "FontSize", 11, "FontWeight", "normal");
            set(ax, "FontSize", 8);

            if row < 3
                set(ax, "XTickLabel", []);
            else
                xlabel("x- axis (cm)");
            end
            if col > 1
                set(ax, "YTickLabel", []);
            else
                ylabel("y- axis (cm)");
            end
            if row == 1 && col == 1
                legend("reference", "simulation", "Location", "northwest", "FontSize", 7);
            end

            validation.patterns(row,col).L = L;
            validation.patterns(row,col).reference = refPanel;
            validation.patterns(row,col).simulation = simPanel;
            validation.patterns(row,col).rmse = rmse;
            validation.patterns(row,col).span = out.span;
        end
    end

    title(tl, "Reference-vs-simulation endpoint trajectories for nine gait envelopes", ...
        "FontSize", 12, "FontWeight", "bold");
end

function [simPanel, refPanel] = panelNormalizeCurves(simCurve, refCurve)
%PANELNORMALIZECURVES Translate curves to a common 0-to-70 cm plotting frame.
    simPanel = simCurve;
    refPanel = refCurve;

    xmin = min([simPanel(1,:), refPanel(1,:)]);
    ymin = min([simPanel(2,:), refPanel(2,:)]);

    simPanel(1,:) = simPanel(1,:) - xmin + 8;
    refPanel(1,:) = refPanel(1,:) - xmin + 8;
    simPanel(2,:) = simPanel(2,:) - ymin + 1.5;
    refPanel(2,:) = refPanel(2,:) - ymin + 1.5;
end

function animateJansenMechanism(main, L, frameStep)
%ANIMATEJANSENMECHANISM Visible one-cycle mechanism simulation.
    if nargin < 3
        frameStep = 3;
    end

    if ~usejava("desktop")
        disp("Live animation skipped because MATLAB is running without the desktop UI.");
        return;
    end

    fig = figure("Name", "Running simulation: modified Jansen mechanism", ...
        "Color", "w", "Position", [160 120 900 560]);
    ax = axes("Parent", fig);

    allX = [main.P0(1,:), main.P1(1,:), main.P2(1,:), main.P3(1,:), main.P4(1,:), ...
        main.P5(1,:), main.P6(1,:), main.PE(1,:)];
    allY = [main.P0(2,:), main.P1(2,:), main.P2(2,:), main.P3(2,:), main.P4(2,:), ...
        main.P5(2,:), main.P6(2,:), main.PE(2,:)];
    xLim = [min(allX)-8, max(allX)+8];
    yLim = [min(allY)-8, max(allY)+8];

    for k = 1:frameStep:numel(main.gaitCycle)
        if ~isvalid(fig)
            break;
        end
        if ~isvalid(ax)
            ax = axes("Parent", fig);
        end
        cla(ax);
        axes(ax); %#ok<LAXES>
        drawMechanism(main, k, [0 0.2 0.85]);
        hold on;
        plot(main.PE(1,1:k), main.PE(2,1:k), "g-", "LineWidth", 2.2);
        plot(main.PE(1,k), main.PE(2,k), "bd", "MarkerFaceColor", "r", "MarkerSize", 7);
        grid on; axis equal;
        xlim(xLim); ylim(yLim);
        xlabel("x position (cm)");
        ylabel("y position (cm)");
        title(sprintf("Running one-DOF simulation, crank angle = %.1f deg", ...
            main.crank(k)*180/pi));
        text(xLim(1)+2, yLim(2)-4, ...
            sprintf("L1=%.2f cm, L4=%.2f cm, L8=%.2f cm", L.val(1), L.val(4), L.val(8)), ...
            "BackgroundColor", "w", "Margin", 4);
        drawnow;
        pause(0.005);
    end
end

function plotEndEffector(main, ref, targetSpan, L)
    figure("Name", "End-effector gait trajectory", "Color", "w");
    plot(main.PE(1,:), main.PE(2,:), "k-", "LineWidth", 2.4); hold on;
    plot(ref(1,:), ref(2,:), "b", "LineWidth", 1.7);
    plot(main.PE(1,1), main.PE(2,1), "ko", "MarkerFaceColor", "k", "MarkerSize", 5);
    grid on; axis equal;
    xlabel("x position (cm)");
    ylabel("y position (cm)");
    title("Modified Jansen end-effector / ankle trajectory");
    legend("simulated mechanism", ...
        sprintf("gait-like reference %.2f x %.2f cm", targetSpan(1), targetSpan(2)), ...
        "start", "Location", "best");

    text(min(main.PE(1,:)), max(main.PE(2,:)), ...
        sprintf("L1=%.2f, L4=%.2f, L8=%.2f cm", L.val(1), L.val(4), L.val(8)), ...
        "VerticalAlignment", "top", "BackgroundColor", "w", "Margin", 4);
end

function plotGaitCurves(main, ref)
    figure("Name", "Gait-cycle position curves", "Color", "w");
    tiledlayout(2,1, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    plot(main.gaitCycle, main.PE(1,:), "r-", "LineWidth", 2.2); hold on;
    plot(main.gaitCycle, ref(1,:), "b--", "LineWidth", 1.5);
    grid on;
    ylabel("x (cm)");
    title("Horizontal ankle motion over gait cycle");
    legend("simulated", "reference envelope", "Location", "best");

    nexttile;
    plot(main.gaitCycle, main.PE(2,:), "c:", "LineWidth", 2.2); hold on;
    plot(main.gaitCycle, ref(2,:), "r--", "LineWidth", 1.5);
    grid on;
    xlabel("gait cycle (%)");
    ylabel("y (cm)");
    title("Vertical ankle motion over gait cycle");
end

function plotVelocityCurves(main)
    figure("Name", "End-effector velocity curves", "Color", "w");
    tiledlayout(3,1, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    plot(main.gaitCycle, main.vx, "LineWidth", 1.8);
    grid on; ylabel("vx (cm/s)");
    title("End-effector velocity with constant crank speed");

    nexttile;
    plot(main.gaitCycle, main.vy, "LineWidth", 1.8);
    grid on; ylabel("vy (cm/s)");

    nexttile;
    plot(main.gaitCycle, main.speed, "LineWidth", 1.8);
    grid on; xlabel("gait cycle (%)"); ylabel("|v| (cm/s)");
end

function plotConstraintResiduals(main)
    figure("Name", "Closed-loop validation", "Color", "w");
    semilogy(main.gaitCycle, main.resnorm + eps, "k-", "LineWidth", 1.8);
    grid on;
    xlabel("gait cycle (%)");
    ylabel("||loop residual||_2 (cm)");
    title("Vector-loop closure residual from fsolve");
end

function plotMechanismSnapshots(main, L)
    figure("Name", "Mechanism snapshots", "Color", "w");
    idx = round(linspace(1, numel(main.gaitCycle)-1, 7));
    for k = idx
        drawMechanism(main, k, [0.65 0.65 0.65]);
        hold on;
    end
    drawMechanism(main, idx(2), [0 0.2 0.8]);
    plot(main.PE(1,:), main.PE(2,:), "m-", "LineWidth", 2.0);
    grid on; axis equal;
    xlabel("x (cm)"); ylabel("y (cm)");
    title(sprintf("Assembly snapshots, L = [%.2f %.1f %.1f %.2f ...] cm", ...
        L.val(1), L.val(2), L.val(3), L.val(4)));
end

function drawMechanism(main, k, color)
    linePairs = {
        "P0","P1"; "P1","P2"; "P3","P2"; "P0","P3";
        "P1","P5"; "P3","P5"; "P2","P4"; "P3","P4";
        "P5","P6"; "P4","P6"; "P5","PE"; "PE","P6"};

    for i = 1:size(linePairs,1)
        A = main.(linePairs{i,1})(:,k);
        B = main.(linePairs{i,2})(:,k);
        plot([A(1), B(1)], [A(2), B(2)], "-", "Color", color, "LineWidth", 2);
        hold on;
    end
    pts = [main.P0(:,k), main.P1(:,k), main.P2(:,k), main.P3(:,k), ...
        main.P4(:,k), main.P5(:,k), main.P6(:,k), main.PE(:,k)];
    plot(pts(1,:), pts(2,:), "d", "Color", color, "MarkerFaceColor", color, "MarkerSize", 5);
end

function sensitivity = runAdjustableLinkSensitivity(Lnom, sim)
%RUNADJUSTABLELINKSENSITIVITY Quantify effects of L1, L4, and L8.
    adjustable = [1 4 8];
    pct = linspace(-0.04, 0.04, 9);
    spans = nan(numel(adjustable), numel(pct), 2);
    areas = nan(numel(adjustable), numel(pct));

    for i = 1:numel(adjustable)
        for j = 1:numel(pct)
            L = Lnom;
            L(adjustable(i)) = Lnom(adjustable(i)) * (1 + pct(j));
            out = simulateJansenCycle(L, sim.NAdvanced, sim.omega, sim.phaseOffset, [], ...
                sim.gaitFrameRotationDeg*pi/180);
            spans(i,j,:) = out.span;
            areas(i,j) = out.area;
        end
    end

    sensitivity = struct();
    sensitivity.adjustable = adjustable;
    sensitivity.percentChange = 100*pct;
    sensitivity.spans = spans;
    sensitivity.areas = areas;
end

function plotSensitivity(sensitivity)
    figure("Name", "Adjustable-link sensitivity", "Color", "w");
    tiledlayout(1,2, "Padding", "compact", "TileSpacing", "compact");
    names = ["L1", "L4", "L8"];

    nexttile;
    for i = 1:3
        plot(sensitivity.percentChange, squeeze(sensitivity.spans(i,:,1)), ...
            "o-", "LineWidth", 1.6); hold on;
    end
    grid on; xlabel("link length change (%)"); ylabel("x-span / stride (cm)");
    title("Effect on horizontal span");
    legend(names, "Location", "best");

    nexttile;
    for i = 1:3
        plot(sensitivity.percentChange, squeeze(sensitivity.spans(i,:,2)), ...
            "o-", "LineWidth", 1.6); hold on;
    end
    grid on; xlabel("link length change (%)"); ylabel("y-span / step height (cm)");
    title("Effect on vertical span");
    legend(names, "Location", "best");
end

function lsDemo = runLeastSquaresSpanMap(Lnom, sim, targetSpan)
%RUNLEASTSQUARESSPANMAP Demonstrates Lambda = Psi*Sigma from the paper.
% Lambda stores [L1; L4; L8] samples. Sigma stores [xspan; yspan] samples.
    L1set = Lnom(1) + [-0.6 0 0.6];
    L4set = Lnom(4) + [-1.0 0 1.0];
    L8set = Lnom(8) + [-1.0 0 1.0];

    n = numel(L1set) * numel(L4set) * numel(L8set);
    Lambda = nan(3, n);
    Sigma = nan(2, n);
    samples = nan(n, 5);
    c = 0;

    for a = 1:numel(L1set)
        for b = 1:numel(L4set)
            for d = 1:numel(L8set)
                c = c + 1;
                L = Lnom;
                L(1) = L1set(a);
                L(4) = L4set(b);
                L(8) = L8set(d);
                out = simulateJansenCycle(L, sim.NAdvanced, sim.omega, sim.phaseOffset, [], ...
                    sim.gaitFrameRotationDeg*pi/180);
                Lambda(:,c) = [L(1); L(4); L(8)];
                Sigma(:,c) = out.span;
                samples(c,:) = [L(1), L(4), L(8), out.span(1), out.span(2)];
            end
        end
    end

    Psi = Lambda * pinv(Sigma);
    LpredAdj = Psi * targetSpan;
    Lpred = Lnom;
    Lpred([1 4 8]) = LpredAdj;
    predicted = simulateJansenCycle(Lpred, sim.N, sim.omega, sim.phaseOffset, [], ...
        sim.gaitFrameRotationDeg*pi/180);

    lsDemo = struct();
    lsDemo.Lambda = Lambda;
    lsDemo.Sigma = Sigma;
    lsDemo.Psi = Psi;
    lsDemo.samples = samples;
    lsDemo.Lpred = Lpred;
    lsDemo.predicted = predicted;
end

function plotLeastSquaresDemo(lsDemo, targetSpan)
    figure("Name", "Least-squares span-to-link map", "Color", "w");
    tiledlayout(1,2, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    scatter(lsDemo.Sigma(1,:), lsDemo.Sigma(2,:), 45, "k", "filled"); hold on;
    plot(targetSpan(1), targetSpan(2), "rp", "MarkerSize", 14, "MarkerFaceColor", "r");
    plot(lsDemo.predicted.span(1), lsDemo.predicted.span(2), "bo", "MarkerSize", 9, "LineWidth", 2);
    grid on;
    xlabel("x-span (cm)");
    ylabel("y-span (cm)");
    title("Sampled span cloud and requested span");
    legend("simulation samples", "requested span", "span from predicted links", "Location", "best");

    nexttile;
    plot(lsDemo.predicted.PE(1,:), lsDemo.predicted.PE(2,:), "b-", "LineWidth", 2.2);
    grid on; axis equal;
    xlabel("x (cm)"); ylabel("y (cm)");
    title(sprintf("LS predicted L1=%.2f, L4=%.2f, L8=%.2f cm", ...
        lsDemo.Lpred(1), lsDemo.Lpred(4), lsDemo.Lpred(8)));
end

function s = spanOf(v)
%SPANOF Max-min span without relying on toolbox-specific range behavior.
    s = max(v(:)) - min(v(:));
end
