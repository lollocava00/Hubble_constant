MODULE funzione
	CONTAINS
	!la funzione da integrare rispetto al redshift
	REAL*8 FUNCTION f(x)
		IMPLICIT NONE
		REAL*8::x
		REAL*8, PARAMETER:: OMEGA_M=0.3d0,OMEGA_LAMBDA=0.7d0
		f=1.d0/SQRT(OMEGA_M*(1+x)**3+OMEGA_LAMBDA)	!<--questa è la funzione che vogliamo integrare
	END FUNCTION f
END MODULE funzione

SUBROUTINE integrale_newton_cotes(integrale,b) 
	USE funzione
	IMPLICIT NONE
	REAL*8::a,delta,x1,x2,summa1,summa2
	INTEGER::n,i
	REAL*8, INTENT(in):: b
	REAL*8, INTENT(out):: integrale

	a=0
	n=10**2			!numero di intervallini
	delta=(b-a)/n 	!spessore intervallini
	integrale=0
	summa1=0
	summa2=0
	!integrale secondo simpson 1/3(ordine 2)
	DO i=1,n
		x1=a+i*delta
		x2=a+(i-1)*delta
		summa1=summa1+f(x1)
		summa2=summa2+f(x2)
	END DO
	integrale=(delta*(4*summa1+2*summa2))/6
END SUBROUTINE integrale_newton_cotes

SUBROUTINE stima_costante_hubble(k,x_int,b_stima,Rv,E_bv,z,step,H0,MJD_max,b_max,b_15,m_15,mu,d_l)
	IMPLICIT NONE
	REAL*8, INTENT(IN OUT):: H0,b_max,MJD_max,m_15,b_15,mu,d_l
	REAL*8:: MJD_15,b_abs,c,b_max_corr,integrale
	INTEGER::j
	INTEGER, INTENT(IN):: k
	REAL*8, DIMENSION(k), INTENT(IN):: b_stima,x_int
	REAL*8, INTENT(IN):: step,Rv,E_bv,z

	b_max=MINVAL(b_stima)
	DO j=1,k
		IF (b_stima(j)== b_max) THEN
			MJD_max=x_int(j)
		END IF
	END DO

	MJD_15=MJD_max+15
	DO j=1,k
		IF(x_int(j)<MJD_15+(step/2) .AND. x_int(j)>MJD_15-(step/2)) THEN
			b_15=b_stima(j)
			MJD_15=x_int(j)
		END IF
	END DO
	!la magnitudine del massimo va corretta per l'estinzione
	b_max_corr=b_max-(Rv+1)*E_bv
	!la differnza tra la mag massima e quella 15 giorni dopo
	m_15=b_15-b_max_corr
	!formula di Hamuy per stimare la magnitudine assoluta al massimo
	b_abs=-19.258d0+0.784d0*(m_15-1.1d0)
	!mu è il modulo di distanza
	mu=b_max_corr-b_abs
	!stima della distanza di luminosità in Mpc
	d_l=(10**((mu+5)/5))/1.0d6

	CALL integrale_newton_cotes(integrale,z)
	c=299792.458d0 !in KM/s
	!stima della costante di hubble
	H0=(c*(1+z)/d_l)*integrale
END SUBROUTINE stima_costante_hubble

SUBROUTINE montecarlo(i,h,j,mean,sigma,z1)
	IMPLICIT NONE
	REAL*8, INTENT(IN):: mean, sigma
	INTEGER, INTENT(IN):: i,h,j !per la determinazione del seed
	INTEGER, PARAMETER:: m=259200, c=54773, a=7141
	INTEGER:: iseed, jran
	REAL*8:: pi,x1,x2,y1,y2
	REAL*8, INTENT(OUT):: z1
	pi=3.14159d0
	iseed=i+h+j
	jran=iseed
	jran=MOD(jran*a+c,m)
	x1=FLOAT(jran)/FLOAT(m)
	jran=MOD(jran*a+c,m)
	x2=FLOAT(jran)/FLOAT(m)
	y1=SQRT(-2.*LOG(x1))*COS(2.*pi*x2)
	!y2=SQRT(-2.*LOG(x1))*SIN(2.*pi*x2)
	z1=sigma*y1+mean
END SUBROUTINE montecarlo

SUBROUTINE spline_cubica(nodi,x,f,stima,k,x_int,step)
	IMPLICIT NONE
	REAL*8::T1,T2,T3,T4,n_max,n_min,fakt,summa
	REAL*8, DIMENSION(nodi), INTENT(in):: x,f
	REAL*8, ALLOCATABLE::d(:),d2_f(:),c(:),a(:),b(:),d2_f2(:),b2(:),d2(:)
	REAL*8, DIMENSION(k), INTENT(out)::stima,x_int
	REAL*8, INTENT(out):: step
	INTEGER, INTENT(in)::nodi,k
	INTEGER::i,n,j,m
	
!la variabile nodi è il numero di punti
	n=nodi-2
	ALLOCATE(d2_f(n),c(n),a(n),b(n),d(n),b2(n),d2(n),d2_f2(nodi))
!la funzione matrice tridiagonale va a creare la matrice nxn le cui soluzioni sono le derivate seconde dei punti x
!diagonallizzando la matrice tridiagonale tramite thomas troviamo i valori delle derivate seconde
!costruzione dei vettori diagonali della matrice (a,b,c)
	DO i=2,n+1
		c(i-1)=x(i)-x(i-1)
		a(i-1)=2*(x(i+1)-x(i-1))
		b(i-1)=x(i+1)-x(i)
		d(i-1)=((6*(f(i+1)-f(i)))/(x(i+1)-x(i)))+((6*(f(i-1)-f(i)))/(x(i)-x(i-1)))
	END DO
!risoluzione matrice tramite algoritmo di thomas
	b2(1)=b(1)/a(1)
	d2(1)=d(1)/a(1)
	DO i=2,n-1
		b2(i)=b(i)/(a(i)-c(i)*b2(i-1))
	END DO
	DO i=2,n
		d2(i)=(d(i)-c(i)*d2(i-1))/(a(i)-c(i)*b2(i-1))
	END DO
	d2_f(n)=d2(n)
	DO i=n-1,1,-1
		d2_f(i)=d2(i)-b2(i)*d2_f(i+1)
	END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!definizione del vettore derivate seconde in cui il primo e ultimo valore sono =0
	d2_f2(1)=0
	d2_f2(nodi)=0
	DO i=1,n
		d2_f2(i+1)=d2_f(i)
	END DO
	n_max = MAXVAL(x)
	n_min = MINVAL(x)
	step = (n_max-n_min)/k
	x_int(1)= n_min
	DO i = 2,k
	   x_int(i) = x_int(i-1) + step
	END DO
!formula dell'interpolazione3
	DO j=1,k
		DO i=1,n+1
			IF(x_int(j)>=x(i).AND.x_int(j)<=x(i+1)) THEN
				T1=(d2_f2(i)*(x(i+1)-x_int(j))**3)/(6*(x(i+1)-x(i)))
				T2=(d2_f2(i+1)*(x_int(j)-x(i))**3)/(6*(x(i+1)-x(i)))
				T3=((f(i)/(x(i+1)-x(i)))-((d2_f2(i)*(x(i+1)-x(i)))/6))*(x(i+1)-x_int(j))
				T4=((f(i+1)/(x(i+1)-x(i)))-((d2_f2(i+1)*(x(i+1)-x(i)))/6))*(x_int(j)-x(i))
				stima(j)=T1+T2+T3+T4
				EXIT
			END IF
		END DO
	END DO
END SUBROUTINE spline_cubica

PROGRAM PROGETTO
	IMPLICIT NONE
	REAL*8::b_max,check_b, MJD_15,MJD_max,b_15,step,m_15,b_max_corr,b_abs,d_l,integrale,c,H0, &
	mu,z1,sum,SIGMA,Hubble_constant,ERR_Hubble,big_sum1,big_sum2,percento, &
	MJD_max2,b_max2,b_152,m_152,mu2,d_l2,SIGMA_MJD_max,SIGMA_b_max,SIGMA_b_15,SIGMA_m_15,SIGMA_mu,SIGMA_d_l, &
	sum2,sum3,sum4,sum5,sum6,sum7
	CHARACTER(6), ALLOCATABLE:: nomi(:)
	REAL*8, ALLOCATABLE::b(:),err_b(:),z(:),nan1(:),nan2(:), b_montecarlo(:),&
	E_bv(:),err_bv(:),Rv(:),err_Rv(:),MJD(:),b2(:),err_b2(:),MJD2(:),b_stima(:),x_int(:), &
	b_stima2(:),H0_montecarlo(:)
	INTEGER::conta,n,i,j,k,h,conta2,conta3,k2,j2,k3,n_file,g,s
	
!ESTRAZIONE DEI DATI DAL FILE newtable_2.txt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	OPEN(1000,file="newtable2.txt")
	conta=0
	DO 
		read(1000,*,END=200)
		conta=conta+1
	END DO
200 continue 
	n=conta-31
	print*, "il file 'newtable2.txt' contiene i dati di: ", n, "supernovae"
	ALLOCATE(nomi(n),z(n),E_bv(n),err_bv(n),Rv(n),err_Rv(n))
	REWIND(1000)
	DO i=1,31
		READ(1000,*) 
	END DO
	DO i=1,n
		READ(1000,"(a6,9x,f7.5,38x,f5.3,x,f5.3,x,f3.1,x,f3.1)") nomi(i), z(i), E_bv(i), err_bv(i), Rv(i), err_Rv(i)
	END DO
	CLOSE(1000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!CORPO DEL PROGRAMMA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	OPEN(1,file='tabella_quantità.dat')
    OPEN(2,file='tabella_H0.dat')
	WRITE(1,'(12x,a7,20x,a5,20x,a4,19x,a10,14x,a2,21x,a20)') 'MJD_max','b_max','b_15','delta_m15','mu','luminosity_distance'
    WRITE(2,'(a4,9x,a2,14x,a5)') 'name', 'H0', 'error'
	PRINT*, "trascrizione dei dati delle quantità: b_max, MJD_max, b_15, m_15, mu, d_l"
	PRINT*, ''
!QUESTO È IL CICLO DO CHE SVOLGERÀ TUTTO IL NECESSARIO PER OGNI SUPERNOVA
	big_sum1=0
	big_sum2=0
	DO i=1,n	
        !ESTRAZIONE DEI DATI DELLE SINGOLE SNe
		n_file=i*10
	 	OPEN(n_file,file='SN'//nomi(i)//'.dat',action='read')
		conta2=0
		DO 
			read(n_file,*,END=300)
			conta2=conta2+1
		END DO
        300 CONTINUE
		k=conta2-5
	    ALLOCATE(MJD(k),nan1(k),nan2(k),b(k),err_b(k))
		REWIND(n_file)
		DO j=1,5
			READ(n_file,*)
		END DO
		!CREAZIONE DI VETTORI CONTENENTI I DATI
		conta3=0
		DO j=1,k
			!nan1 e nan2 sono quantità che non sono utili allo svolgimento del progetto
			READ(n_file,*) MJD(j),nan1(j),nan2(j), b(j), err_b(j)
			IF (b(j)>90.0d0) THEN
				conta3=conta3+1
			END IF
		END DO
		CLOSE(n_file)
		!PULIZIA DELLE MAG CON DATI 99.900 TRAMITE CREAZIONE DI NUOVI VETTORI: b2, MJD2, err_b2
		k2=k-conta3
		ALLOCATE(b2(k2),MJD2(k2),err_b2(k2))
		j2=1
		DO j=1,k
			IF (b(j)<90.0d0) THEN
				b2(j2)=b(j)
				MJD2(j2)=MJD(j)
				err_b2(j2)=err_b(j)
				j2=j2+1
			END IF
		END DO

		k3= 50*size(MJD2)	!!!!!!!!!!!!!!!!!!!!!!!!!!!!numero di valori di MJD su cui voglio interpolare

		ALLOCATE(b_stima(k3),x_int(k3))
		!CONTROLLO CHE NON CI SIANO SOLO CURVE DI LUCE MONOTONE
		check_b=MINVAL(b2)
		IF (check_b==b2(1)) THEN
			GO TO 400
		END IF
		!INTERPOLAZIONE CON SPLINE PER OTTENERE UNA STIMA DELLA CURVA DI LUCE NEL MASSIMO DI MAG.
		CALL spline_cubica(k2,MJD2,b2,b_stima,k3,x_int,step)
		!STIMA DELLA COSTANTE DI HUBBLE(CHE SARA FATTA PER OGNI SUPERNOVA)
		CALL stima_costante_hubble(k3,x_int,b_stima,Rv(i),E_bv(i),z(i),step,H0, &
		MJD_max,b_max,b_15,m_15,mu,d_l)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!CALCOLO DEGLI ERRORI TRAMITE GENERAZIONE DI NUMERI CASUALI DISTRIBUITI GAUSSIANAMENTE!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		ALLOCATE(H0_montecarlo(100))
		sum=0.d0
		sum2=0.d0
		sum3=0.d0
		sum4=0.d0
		sum5=0.d0
		sum6=0.d0
		sum7=0.d0
		DO h=1,100
			ALLOCATE(b_montecarlo(k2),b_stima2(k3))
			DO j=1,k2
				CALL montecarlo(i,h,j,b2(j),err_b2(j),z1)
				b_montecarlo(j)=z1
				!print*, b_montecarlo(j)
			END DO
			CALL spline_cubica(k2,MJD2,b_montecarlo,b_stima2,k3,x_int,step)
			CALL stima_costante_hubble(k3,x_int,b_stima2,Rv(i),E_bv(i),z(i),step,H0_montecarlo(h), &
			MJD_max2,b_max2,b_152,m_152,mu2,d_l2)
            !queste somme vanno a calcolare gli errori per ogni quantità richiesta
			sum=sum+(H0_montecarlo(h)-H0)**2
			sum2=sum2+(MJD_max2-MJD_max)**2
			sum3=sum3+(b_max2-b_max)**2
			sum4=sum4+(b_152-b_15)**2
			sum5=sum5+(m_152-m_15)**2
			sum6=sum6+(mu2-mu)**2
			sum7=sum7+(d_l2-d_l)**2
			
			DEALLOCATE(b_montecarlo,b_stima2)
		END DO
		DEALLOCATE(H0_montecarlo)
		SIGMA=SQRT(sum/100)
		! big_sum SONO LE SOMME PER LA MEDIA PESATA
		big_sum1=big_sum1+H0/(SIGMA**2)
		big_sum2=big_sum2+1/(SIGMA**2)
		PRINT*, 'SN',nomi(i),' H0:', H0, '+/-', SIGMA, '(Km/s/Mpc)'
		!CREAZIONE DATI DEGLI ERRORI DA METTERE NELLE TABELLE
		SIGMA_MJD_max=SQRT(sum2/100)
		SIGMA_b_max=SQRT(sum3/100)
		SIGMA_b_15=SQRT(sum4/100)
		SIGMA_m_15=SQRT(sum5/100)
		SIGMA_mu=SQRT(sum6/100)
		SIGMA_d_l=SQRT(sum7/100)
		WRITE(1,'(a2,a6,4x,f10.4,2x,a3,2x,f7.4,4x,f7.4,2x,a3,2x,f6.4, &
		4x,f7.4,2x,a3,2x,f6.4,4x,f6.4,2x,a3,2x,f6.4, &
		4x,f7.4,2x,a3,2x,f6.4,4x,f8.4,2x,a3,2x,f7.4)') 'SN',nomi(i) ,mjd_max, '+/-', &
		SIGMA_mjd_max, b_max, '+/-', SIGMA_b_max, b_15, '+/-', SIGMA_b_15, m_15, '+/-', SIGMA_m_15, &
		mu, '+/-', SIGMA_mu, d_l, '+/-', SIGMA_d_l
        WRITE(2,'(a2,a6,4x,f10.6,2x,a3,2x,f8.6)') 'SN',nomi(i), H0, '+/-',SIGMA
400 CONTINUE
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DEALLOCATE(MJD,nan1,nan2,b,err_b,MJD2,b2,err_b2,b_stima,x_int)
	END DO

	CLOSE(1)
    CLOSE(2)

	PRINT*, ''
    PRINT*, "i risultati sono stati salvati nel file: 'tabella_H0.dat' "
	PRINT*, "i risultati sono stati salvati nel file: 'tabella_quantità.dat'"
	PRINT*, ''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	PRINT*, 'Estrazione dei dati per ogni supernova Completato.'
	PRINT*, 'Il programma ha svolto il calcolo della costante di hubble di tutte le supernove.'
	!ORA CHE ABBIAMO UN VALORE DI H0 E IL SUO ERRORE PER OGNI SUPERNOVA
	!FACCIAMO UNA MEDIA PESATA PER OTTENERE UNA STIMA FINALE DELLA COSTANTE DI HUBBLE
	Hubble_constant=big_sum1/big_sum2
	ERR_Hubble=SQRT(1/big_sum2)
	PRINT*, 'la stima finale della costante di hubble è:', Hubble_constant, '+/-', ERR_Hubble, '(Km/s/Mpc)'
END PROGRAM PROGETTO