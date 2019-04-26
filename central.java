
import java.io.*;
import java.util.*;

import org.jlab.groot.math.*;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.base.GStyle;

public class central {
	boolean userTimeBased, write_volatile;
	public int runNum;
	public boolean BackToBack;
	public float STT, RFT, MinCTOF,MaxCTOF, minSTT, maxSTT;
	public float[] CTOF_shft;

	public H2F H_CTOF_pos, H_CTOF_edep_phi, H_CTOF_edep_z, H_CTOF_path_mom;
	public H2F H_CTOF_edep_pad_neg, H_CTOF_edep_pad_pos;
	public H2F H_vz_DC_CVT, H_phi_DC_CVT, H_CVT_CTOF_phi, H_CVT_CTOF_z, H_CVT_t_STT, H_CVT_t_pad;
	public H1F[] H_CVT_t;
	public H1F H_CVT_t_pos, H_CVT_t_neg;

	//for timeline
	public H1F H_CTOF_pos_mass, H_CTOF_neg_mass;
	public H2F H_CTOF_vt_pim;
	public H1F H_CTOF_edep_pim;

	public central(int reqrunNum, boolean reqTimeBased, boolean reqwrite_volatile) {
		runNum = reqrunNum;userTimeBased=reqTimeBased;
		write_volatile = reqwrite_volatile;
                CTOF_shft = new float[]{0.00f , -300.51f , -297.72f , -295.05f , -301.03f , -298.19f , -295.41f , -299.74f , -298.30f , -295.21f ,
                                        -300.22f , -297.35f , -295.28f , -299.14f , 0.00f , -294.51f , -298.80f , -298.30f , -295.49f , -298.81f ,
                                        -297.23f , -294.44f , -298.96f , -296.27f , -295.30f , -298.74f , -296.44f , -293.98f , -299.07f , -296.08f ,
                                        -293.52f , -297.96f , -296.42f , -294.40f , -298.90f , -297.14f , -295.28f , -299.03f , -297.71f , -294.35f ,
                                        -298.81f , -296.75f , -293.77f , -298.96f , -297.65f , -294.67f , -299.73f , -297.76f , -295.87f , -297.10f};//2476

		H_CTOF_pos = new H2F("H_CTOF_pos","H_CTOF_pos",50,-180,180,100,-10,5);
		//H_CTOF_pos = new H2F("H_CTOF_pos","H_CTOF_pos",50,-180,180,100,-20,10);
		H_CTOF_pos.setTitle("CTOF hits z vs #phi");
		H_CTOF_pos.setTitleX("#phi (^o)");
		H_CTOF_pos.setTitleY("z/10 (cm)");
		H_CTOF_edep_phi = new H2F("H_CTOF_edep_phi","H_CTOF_edep_phi",50,-180,180,100,0,60);
		H_CTOF_edep_phi.setTitle("CTOF Edep vs #phi");
		H_CTOF_edep_phi.setTitleX("#phi (^o)");
		H_CTOF_edep_phi.setTitleY("E (MeV)");
		H_CTOF_edep_pad_pos = new H2F("H_CTOF_edep_phi","H_CTOF_edep_phi",50,0.5,50.5,200,0,60);
		H_CTOF_edep_pad_pos.setTitle("CTOF Edep vs pad, pos. tracks");
		H_CTOF_edep_pad_pos.setTitleX("pad");
		H_CTOF_edep_pad_pos.setTitleY("E (MeV)");
                H_CTOF_edep_pad_neg = new H2F("H_CTOF_edep_phi","H_CTOF_edep_phi",50,0.5,50.5,200,0,60);
                H_CTOF_edep_pad_neg.setTitle("CTOF Edep vs pad, neg. tracks");
                H_CTOF_edep_pad_neg.setTitleX("pad");
                H_CTOF_edep_pad_neg.setTitleY("E (MeV)");
		H_CTOF_edep_z = new H2F("H_CTOF_edep_z","H_CTOF_edep_z",100,-10,5,100,0,150);
		H_CTOF_edep_z.setTitle("CTOF Edep vs z");
		H_CTOF_edep_z.setTitleX("z/10 (cm)");
		H_CTOF_edep_z.setTitleY("E (MeV)");
		H_vz_DC_CVT = new H2F("H_vz_DC_CVT","H_vz_DC_CVT",100,-20,20,100,-20,20);
		H_vz_DC_CVT.setTitle("CVT vz vs DC vz");
		H_vz_DC_CVT.setTitleX("DC vz (cm)");
		H_vz_DC_CVT.setTitleY("CVT vz (cm)");
		H_phi_DC_CVT = new H2F("H_phi_DC_CVT","H_phi_DC_CVT",100,-180,180,100,-180,180);
		H_phi_DC_CVT.setTitle("CVZ #phi vs DC #phi");
		H_phi_DC_CVT.setTitleX("DC #phi (^o)");
		H_phi_DC_CVT.setTitleY("CVT #phi (^o)");
		H_CTOF_path_mom = new H2F("H_CTOF_path_mom","H_CTOF_path_mom",100,0,2,100,20,60);
		H_CTOF_path_mom.setTitle("CTOF path vs mom");
		H_CTOF_path_mom.setTitleX("p (GeV)");
		H_CTOF_path_mom.setTitleY("CTOF path (cm)");
		H_CVT_CTOF_phi = new H2F("H_CVT_CTOF_phi","H_CVT_CTOF_phi",100,-180,180,50,-180,180);
		H_CVT_CTOF_phi.setTitle("CVT CTOF #Delta #phi");
		H_CVT_CTOF_phi.setTitleX("CVT #phi (^o)");
		H_CVT_CTOF_phi.setTitleY("CTOF #phi (^o)");
		H_CVT_CTOF_z = new H2F("H_CVT_CTOF_z","H_CVT_CTOF_z",100,-10,5,100,-10,5);
		H_CVT_CTOF_z.setTitle("CVT CTOF #Delta z");
		H_CVT_CTOF_z.setTitleX("CVT z/10 (cm)");
		H_CVT_CTOF_z.setTitleY("CTOF z/10 (cm)");
		minSTT = 100;maxSTT=200;
		if(runNum>0&&runNum<3210){
			minSTT = 540;maxSTT=590;
		}
		MinCTOF=320;MaxCTOF=400;
		if(runNum>0&&runNum<3210){
			MinCTOF=-5;MaxCTOF=5;
		}
		MinCTOF=-5;MaxCTOF=5;
		H_CVT_t_STT = new H2F("H_CVT_t_STT","H_CVT_t_STT",100,minSTT,maxSTT,100,minSTT+(MinCTOF+MaxCTOF)/2,maxSTT+(MinCTOF+MaxCTOF)/2);
		H_CVT_t_STT.setTitle("CTOF vertex t vs STT, neg. tracks");
		H_CVT_t_STT.setTitleX("STT (ns)");
		H_CVT_t_STT.setTitleY("CTOF vertex t (ns)");
                try {
                        Thread.sleep(5000);// in ms
                }catch (Exception e) {
                        System.out.println(e);
                }
		System.out.println("CTOF range "+MinCTOF+" to "+MaxCTOF);
		H_CVT_t_pad = new H2F("H_CVT_t_pad","H_CVT_t_pad",50,0.5,50.5,250,MinCTOF,MaxCTOF);
		H_CVT_t_pad.setTitle("CTOF vertex t - STT  vs pad, neg. tracks");
		H_CVT_t_pad.setTitleX("pad");
		H_CVT_t_pad.setTitleY("CTOF vertex t - STT (ns)");
		H_CVT_t = new H1F[50];
		for(int p=0;p<50;p++){
			H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),250,MinCTOF,MaxCTOF);
			H_CVT_t[p].setTitle(String.format("pad %d CTOF vertex t - STT, neg. tracks",p+1));
			H_CVT_t[p].setTitleX("Delta_t (ns)");
		}
                H_CVT_t_pos = new H1F("H_CVT_t_pos","H_CVT_t_pos",250,MinCTOF,MaxCTOF);
                H_CVT_t_pos.setTitle("All CTOF pads, CTOF vertex t - STT, pos. tracks");
                H_CVT_t_pos.setTitleX("CTOF vertex t - STT (ns)");

                H_CVT_t_neg = new H1F("H_CVT_t_neg","H_CVT_t_neg",250,MinCTOF,MaxCTOF);
                H_CVT_t_neg.setTitle("All CTOF pads, CTOF vertex t - STT, neg. tracks");
                H_CVT_t_neg.setTitleX("CTOF vertex t - STT (ns)");

		H_CTOF_pos_mass = new H1F("H_CTOF_pos_mass","H_CTOF_pos_mass",100,-0.5,3.5);
		H_CTOF_pos_mass.setTitle("pos Mass^2");
		H_CTOF_pos_mass.setTitleX("M^2 (GeV^2)");
		H_CTOF_neg_mass = new H1F("H_CTOF_neg_mass","H_CTOF_neg_mass",100,-0.5,2.0);
		H_CTOF_neg_mass.setTitle("neg Mass^2");
		H_CTOF_neg_mass.setTitleX("M^2 (GeV^2)");
		H_CTOF_vt_pim = new H2F("H_CTOF_vt_pim","H_CTOF_vt_pim",250,-5,5,250,-5,5);
		H_CTOF_vt_pim.setTitle("CTOF MIP (pi-) vertex time");
		H_CTOF_vt_pim.setTitleX("vertex time - RFTime (ns)");
		H_CTOF_vt_pim.setTitleY("vertex time - STTime (ns)");
		H_CTOF_edep_pim = new H1F("H_CTOF_edep_pim","H_CTOF_edep_pim",100,0,150);
		H_CTOF_edep_pim.setTitle("CTOF MIP (pi-) Edep");
		H_CTOF_edep_pim.setTitleX("E (MeV)");
	}
	public double Vangle(Vector3 v1, Vector3 v2){ 
		double res = 0; 
		double l1 = v1.mag();
		double l2 = v2.mag();
		double prod = v1.dot(v2);
		if( l1 * l2 !=0 && Math.abs(prod)<l1*l2 )res = Math.toDegrees( Math.acos(prod/(l1*l2) ) ); 
		return res; 
	}
	public void FillTracks(DataBank DCbank, DataBank CVTbank){
		for(int iDC=0;iDC<DCbank.rows();iDC++){
			for(int iCVT=0;iCVT<CVTbank.rows() && !BackToBack;iCVT++){
				float DelPhi = (float)Math.toDegrees( Math.atan2(DCbank.getFloat("p0_y",iDC),DCbank.getFloat("p0_x",iDC)) - CVTbank.getFloat("phi0",iCVT));
				//DelPhi += 180f - 10f;
				while(DelPhi>180)DelPhi-=360;
				while(DelPhi<-180)DelPhi+=360;
				if(Math.abs(DelPhi)<180){
					H_vz_DC_CVT.fill(DCbank.getFloat("Vtx0_z",iDC),CVTbank.getFloat("z0",iCVT));
					//H_vz_DC_CVT.fill(DCbank.getFloat("Vtx0_z",iDC),0.1f*CVTbank.getFloat("z0",iCVT));
					H_phi_DC_CVT.fill(Math.toDegrees(Math.atan2(DCbank.getFloat("p0_y",iDC),DCbank.getFloat("p0_x",iDC))),Math.toDegrees(CVTbank.getFloat("phi0",iCVT)));
					BackToBack = true;
				}
			}
		}
	}
	public void FillCVTCTOF(DataBank CVTbank, DataBank CTOFbank){
		for(int iCTOF=0;iCTOF<CTOFbank.rows();iCTOF++){
			float e = CTOFbank.getFloat("energy",iCTOF);
			int trackid = CTOFbank.getInt("trkID",iCTOF);
			if(!Float.isNaN(e) && trackid>-1 && e>1.5){
				boolean matched = false;
				for(int iCVT=0;iCVT<CVTbank.rows() && !matched;iCVT++){
					float mom = CVTbank.getFloat("p",iCVT);
					float cx = CVTbank.getFloat("c_x",iCVT)*0.1f;
					float cy = CVTbank.getFloat("c_y",iCVT)*0.1f;
					float cz = CVTbank.getFloat("c_z",iCVT)*0.1f;
					float cphi = (float)Math.toDegrees(Math.atan2(cy,cx));
					int charge = CVTbank.getInt("q",iCVT);
					int pad = CTOFbank.getInt("component",iCTOF);
					float x = CTOFbank.getFloat("x",iCTOF)*0.1f;
					float y = CTOFbank.getFloat("y",iCTOF)*0.1f;
					float z = CTOFbank.getFloat("z",iCTOF)*0.1f;
					float t = CTOFbank.getFloat("time",iCTOF);
					float p = CTOFbank.getFloat("pathLength",iCTOF);
					float phi = (float)Math.toDegrees(Math.atan2(y,x));
					//float beta =  mom/(float)Math.sqrt(mom*mom+0.93827f*0.93827f);
					float beta =  mom/(float)Math.sqrt(mom*mom+0.13957061f*0.13957061f);
					//float DelPhi = phi-cphi+190;
					float DelPhi = phi-cphi;

					//like FTOF in dst_mon
					double CTOFbeta = p/(29.98f*(t-STT));
					double CTOFmass = mom * mom * ( 1/(CTOFbeta*CTOFbeta) - 1);

					while(DelPhi>180)DelPhi-=360;
					while(DelPhi<-180)DelPhi+=360;
					if( Math.abs(DelPhi)<30 && z<4.7){}
					if( Math.abs(DelPhi)<180){
						H_CVT_CTOF_phi.fill(cphi,phi);
						H_CVT_CTOF_z.fill(cz,z);
						H_CTOF_pos.fill(phi,z);
						H_CTOF_edep_phi.fill(phi,e);
						H_CTOF_edep_z.fill(z,e);
						H_CTOF_path_mom.fill(mom,p);
						float CTOFTime = t - p/29.92f/beta;
						//float CTOFTime = t - p/29.92f/beta - CTOF_shft[pad];
						if (charge < 0) {
							H_CVT_t_STT.fill(STT,CTOFTime);
							H_CVT_t_pad.fill(pad,CTOFTime-STT);
							H_CVT_t[pad].fill(CTOFTime-STT);
							H_CVT_t[49].fill(CTOFTime-STT);
							H_CTOF_edep_pad_neg.fill(pad,e);
							H_CVT_t_neg.fill(CTOFTime-STT);
							H_CTOF_neg_mass.fill(CTOFmass);
							//pi- fiducial cut borrowing from Pierre's CND
							if (Math.sqrt(Math.abs(CTOFmass))<0.38 && CTOFmass>-0.35*0.35){
								H_CTOF_vt_pim.fill(CTOFTime-RFT,CTOFTime-STT);
								H_CTOF_edep_pim.fill(e);
							}			
						}
						if (charge>0) {
							H_CTOF_edep_pad_pos.fill(pad,e);
							H_CVT_t_pos.fill(CTOFTime-STT);
							H_CTOF_pos_mass.fill(CTOFmass);
						}
						matched = true;
					}
				}
			}
		}
	}
	public void processEvent(DataEvent event) {
		BackToBack = false;
		DataBank eventBank = null, trackDetBank = null;
                if(userTimeBased){
                        if(event.hasBank("REC::Event"))eventBank = event.getBank("REC::Event");
                        if(event.hasBank("TimeBasedTrkg::TBTracks"))trackDetBank = event.getBank("TimeBasedTrkg::TBTracks");
                }
                if(!userTimeBased){
                        if(event.hasBank("RECHB::Event"))eventBank = event.getBank("RECHB::Event");
                        if(event.hasBank("HitBasedTrkg::HBTracks"))trackDetBank = event.getBank("HitBasedTrkg::HBTracks");
                }

		if(eventBank!=null)STT = eventBank.getFloat("STTime",0);
		if(eventBank!=null)RFT = eventBank.getFloat("RFTime",0);
		else return;
		if(trackDetBank != null && event.hasBank("CVTRec::Tracks"))FillTracks(trackDetBank,event.getBank("CVTRec::Tracks"));
		if(BackToBack && event.hasBank("CVTRec::Tracks") && event.hasBank("CTOF::hits"))FillCVTCTOF(event.getBank("CVTRec::Tracks"),event.getBank("CTOF::hits"));
	}
	public void plot() {
		EmbeddedCanvas can_central  = new EmbeddedCanvas();
                can_central.setSize(5000,3000);
                can_central.divide(5,3);
                can_central.setAxisTitleSize(30);
                can_central.setAxisFontSize(30);
                can_central.setTitleSize(30);
		can_central.cd(0);can_central.draw(H_CTOF_pos);
		can_central.cd(1);can_central.draw(H_CVT_CTOF_phi);
		can_central.cd(2);can_central.draw(H_CVT_CTOF_z);
		can_central.cd(3);can_central.draw(H_CTOF_path_mom);
		can_central.cd(4);can_central.draw(H_CTOF_edep_phi);
		can_central.cd(5);can_central.draw(H_CTOF_edep_z);
		can_central.cd(6);can_central.draw(H_CTOF_edep_pad_pos);
		can_central.cd(7);can_central.draw(H_CTOF_edep_pad_neg);
		can_central.cd(8);can_central.draw(H_vz_DC_CVT);
		can_central.cd(9);can_central.draw(H_phi_DC_CVT);
		can_central.cd(10);can_central.draw(H_CVT_t_STT);
		can_central.cd(11);can_central.draw(H_CVT_t_pad);

		can_central.cd(12);can_central.draw(H_CVT_t_pos);
                can_central.cd(13);can_central.draw(H_CVT_t_neg);

		if(runNum>0){
			if(!write_volatile)can_central.save(String.format("plots"+runNum+"/central.png"));
			if(write_volatile)can_central.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/central.png"));
			System.out.println(String.format("saved plots"+runNum+"/central.png"));
		}
		else{
			can_central.save(String.format("plots/central.png"));
			System.out.println(String.format("saved plots/central.png"));
		}
		for(int p=0;p<10;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
		for(int p=10;p<20;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
		for(int p=20;p<30;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
		for(int p=30;p<40;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
		for(int p=40;p<50;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
	}
        public void write(){
                TDirectory dirout = new TDirectory();
                dirout.mkdir("/ctof/");
                dirout.cd("/ctof/");
                dirout.addDataSet(H_CVT_t_pad,H_CTOF_edep_phi);
                for(int p=0;p<50;p++)dirout.addDataSet(H_CVT_t[p]);
                dirout.addDataSet(H_CVT_t_pos, H_CVT_t_neg);
                dirout.addDataSet(H_CTOF_pos_mass, H_CTOF_neg_mass, H_CTOF_vt_pim, H_CTOF_edep_pim);

		if(write_volatile)if(runNum>0)dirout.writeFile("/volatile/clas12/rga/spring18/plots"+runNum+"/out_CTOF_"+runNum+".hipo");
                
		if(!write_volatile){
			if(runNum>0)dirout.writeFile("plots"+runNum+"/out_CTOF_"+runNum+".hipo");
			else dirout.writeFile("plots/out_CTOF.hipo");
		}
        }   

////////////////////////////////////////////////
        public static void main(String[] args) {
                System.setProperty("java.awt.headless", "true");
		GStyle.setPalette("kRainBow");
                int count = 0;
		int runNum = 0;
		boolean useTB = true;
		boolean useVolatile = false;
		String filelist = "list_of_files.txt";
		if(args.length>0)runNum=Integer.parseInt(args[0]);
		if(args.length>1)filelist = args[1];
		if(args.length>2)if(Integer.parseInt(args[2])==0)useTB=false;
		central ana = new central(runNum,useTB,useVolatile);
                List<String> toProcessFileNames = new ArrayList<String>();
                File file = new File(filelist);
                Scanner read;
                try {
                        read = new Scanner(file);
                        do { 
                                String filename = read.next();
                                toProcessFileNames.add(filename);

                        }while (read.hasNext());
                        read.close();
                }catch(IOException e){
                        e.printStackTrace();
                }
		//int maxevents = 50000000;
		int maxevents = 500000;
		if(args.length>2)maxevents=Integer.parseInt(args[2]);
		int progresscount=0;int filetot = toProcessFileNames.size();
		for (String runstrg : toProcessFileNames) if( count<maxevents ){
			progresscount++;
			System.out.println(String.format(">>>>>>>>>>>>>>>> %s",runstrg));
			File varTmpDir = new File(runstrg);
			if(!varTmpDir.exists()){System.out.println("FILE DOES NOT EXIST");continue;}
			System.out.println("READING NOW "+runstrg);
			HipoDataSource reader = new HipoDataSource();
			reader.open(runstrg);
			int filecount = 0;
			while(reader.hasEvent()&& count<maxevents ) {
				DataEvent event = reader.getNextEvent();
				ana.processEvent(event);
				filecount++;count++;
				if(count%10000 == 0) System.out.println(count/1000 + "k events (this is central analysis on "+runstrg+") ; progress : "+progresscount+"/"+filetot);
			}
			reader.close();
		}
		System.out.println("Total : " + count + " events");
		ana.plot();
		ana.write();
        }
}
