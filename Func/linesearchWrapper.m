function [actU,actW,sols] = linesearchWrapper(iniActU,iniActW,theta0,thetaE,s)


actU = iniActU;
actW = iniActW;

d = thetaE-theta0;
thetaC = theta0;

iter = 0;
stepsize = 0;

sols = cell(2*s.N^2,1);

fprintf('Iteration : Parameter : Input constraints : Disturbance constraints\n');
fprintf('%2.0d. : %9.2e : [ %8.4f , %8.4f , %8.4f , %8.4f ]  : ', iter,stepsize,thetaC(1),thetaC(2),thetaC(3),thetaC(4))
    for l = 1:2
        for m = 1:s.N
            if l == 1
                fprintf('%d',actU(m,:))
                fprintf(' ')
            elseif l == 2
                fprintf('%d',actW(m,:))
                fprintf(' ')
            end
        end
        fprintf(' : ')
    end
    iter = iter+1;


while d'*thetaC<d'*d-s.tol;

    if iter==6
        0;
    end

[wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,...
    kappaK,kappak,kappaZ,rhoK,rhok,rhoZ] = solver(actU,actW,s);

sols{iter} = struct('wM',wM,'wm',wm,'xpM',xpM,'xpm',xpm,'lambdaM',lambdaM,...
    'lambdam',lambdam,'lambdaZ',lambdaZ,'muM',muM,'mum',mum,'muZ',muZ,...
    'zetaM',zetaM,'zetam',zetam,'zetaZ',zetaZ,'uK',uK,'uk',uk,...
    'kappaK',kappaK,'kappak',kappak,'kappaZ',kappaZ,'rhoK',rhoK,...
    'rhok',rhok,'rhoZ',rhoZ);

[sols{iter}.x,sols{iter}.u,sols{iter}.w,sols{iter}.kappa,sols{iter}.lambda]...
    = trajectoryEvaluator(thetaC,wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,...
    mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,kappaZ,rhoK,rhok,rhoZ,s);

[uIdx,wIdx,stepsize] = linesearch(thetaC,thetaE,actU,actW,wM,wm,xpM,...
    xpm,lambdaM,lambdam,lambdaZ,muM,mum,muZ,zetaM,zetam,zetaZ,uK,uk,...
    kappaK,kappak,kappaZ,rhoK,rhok,rhoZ,s);


    if d'*(thetaC + stepsize*(thetaE-thetaC))>=thetaE'*thetaE-s.tol
        thetaC = thetaE;
        [sols{iter+1}.x,sols{iter+1}.u,sols{iter+1}.w,sols{iter+1}.kappa,sols{iter+1}.lambda]...
        = trajectoryEvaluator(thetaC,wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,...
        mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,kappaZ,rhoK,rhok,rhoZ,s);
        display('Terminated')
        break;
    end

    if and(stepsize<1,stepsize>0)
        actW = wIdx;
        actU = uIdx;
        thetaC = thetaC + stepsize*(thetaE-thetaC);
            fprintf('%2.0d. : %9.2e : [ %8.4f , %8.4f , %8.4f , %8.4f ]  : ', iter,stepsize,thetaC(1),thetaC(2),thetaC(3),thetaC(4))
            for l = 1:2
                for m = 1:s.N
                    if l == 1
                        fprintf('%d',actU(m,:))
                        fprintf(' ')
                    elseif l == 2
                        fprintf('%d',actW(m,:))
                        fprintf(' ')
                    end
                end
                fprintf(' : ')
            end
            iter = iter+1;
    elseif stepsize>=1
        thetaC = thetaE;
        [sols{iter+1}.x,sols{iter+1}.u,sols{iter+1}.w,sols{iter+1}.kappa,sols{iter+1}.lambda]...
        = trajectoryEvaluator(thetaC,wM,wm,xpM,xpm,lambdaM,lambdam,lambdaZ,muM,...
        mum,muZ,zetaM,zetam,zetaZ,uK,uk,kappaK,kappak,kappaZ,rhoK,rhok,rhoZ,s);
        display('Terminated')
        break;
    else
        display('Infeasible problem')
        break;
    end

    
end
fprintf('\n')