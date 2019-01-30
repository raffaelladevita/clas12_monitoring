
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

public class cnd {
	boolean userTimeBased, write_volatile;
	public int runNum;
	public float STT;
	public float[] Tmean;
	
	public H1F H_CND_occ, H_CND_time;
	public H1F[] H_CND_t;
	public H2F H_CND_phi_pad, H_CND_layer_pad;
	public H2F[] H_CND_edep_z, H_CND_edep_phi, H_CND_vt_pad, H_CND_t_edep;
	public H2F[] H_CVT_CND_z, H_CVT_CND_phi;
	//public H2F[] H_neut_Edep_z, H_neut_z_phi, H_neut_bet_pad;
	public cnd(int reqrunNum, boolean reqTimeBased, boolean reqwrite_volatile) {
		userTimeBased=reqTimeBased;
		write_volatile = reqwrite_volatile;
		Tmean = new float[]{
		379.97f , 380.71f , 381.20f , 381.62f , 380.56f , 381.19f , 380.47f , 382.03f , 380.62f , 381.68f , 381.36f , 380.70f , 378.82f , 379.98f , 381.41f , 380.57f , 381.37f , 380.68f , 379.62f , 380.95f ,
		381.58f , 381.02f , 381.44f , 381.17f , 0.00f , 379.96f , 380.14f , 381.30f , 381.11f , 380.59f , 381.11f , 380.83f , 381.48f , 380.50f , 381.03f , 380.73f , 380.37f , 379.71f , 380.08f , 381.91f , 
		379.82f , 382.43f , 381.15f , 380.00f , 380.68f , 381.56f , 381.13f , 381.62f , 381.74f , 0.00f , 379.98f , 380.61f , 381.10f , 381.67f , 381.95f , 380.91f , 380.74f , 379.92f , 381.69f , 381.95f , 
		380.47f , 379.75f , 378.13f , 380.66f , 381.02f , 381.75f , 381.76f , 381.37f , 380.52f , 381.31f , 382.04f , 381.31f , 381.68f , 381.31f , 0.00f , 380.71f , 380.69f , 380.48f , 381.92f , 381.15f , 
		381.48f , 380.95f , 380.06f , 381.17f , 381.65f , 380.25f , 379.82f , 378.51f , 381.16f , 381.34f , 381.29f , 380.94f , 381.52f , 380.75f , 381.82f , 381.66f , 381.02f , 381.26f , 381.31f , 0.00f , 
		382.20f , 381.77f , 381.05f , 381.07f , 381.26f , 381.42f , 381.12f , 381.66f , 381.34f , 381.68f , 380.50f , 379.16f , 380.96f , 381.50f , 380.68f , 380.00f , 381.37f , 380.98f , 380.80f , 380.59f ,
		381.13f , 380.11f , 381.35f , 381.39f , 0.00f , 381.54f , 381.48f , 380.41f , 380.42f , 381.80f , 381.08f , 380.84f , 382.04f , 381.95f , 381.60f , 379.73f , 379.75f , 380.78f , 381.16f , 380.16f , 
		379.51f , 380.56f , 380.19f , 379.90f , 380.46f , 381.52f , 380.25f , 380.54f , 380.91f , 0.00f
		};
		//float T0 = 360;
		//float T1 = 400;
		float T0 = -20;
		float T1 =  20;
		runNum = reqrunNum;
		H_CND_time = new H1F("H_CND_time","H_CND_time",100,T0,T1);
		H_CND_time.setTitle("CND vertex time - STT");
		H_CND_time.setTitleX("vt - STT (ns)");
		H_CND_occ = new H1F("H_CND_occ","H_CND_occ",48,0.5,48.5);
		H_CND_occ.setTitle("CND occupancy");
		H_CND_occ.setTitleX("CND counter #");
		H_CND_phi_pad = new H2F("H_CND_phi_pad","H_CND_phi_pad",48,0.5,48.5,48,-180,180);
		H_CND_phi_pad.setTitle("CND #phi vs pad");
		H_CND_phi_pad.setTitleX("CND counter #");
		H_CND_phi_pad.setTitleY("CND #phi (^o)");
		H_CND_layer_pad = new H2F("H_CND_layer_pad","H_CND_layer_pad",3,0.5,3.5,48,0.5,48.5);
		H_CND_layer_pad.setTitle("CND layer vs pad");
		H_CND_layer_pad.setTitleX("CND counter #");
		H_CND_layer_pad.setTitleY("CND layer");
		H_CND_edep_z = new H2F[3];
		H_CND_edep_phi = new H2F[3];
		H_CVT_CND_z = new H2F[3];
		H_CVT_CND_phi = new H2F[3];
		H_CND_vt_pad  = new H2F[3];
		H_CND_t_edep = new H2F[150];
		H_CND_t = new H1F[150];
		for(int p=0;p<150;p++){
			H_CND_t_edep[p] = new H2F(String.format("H_CND_t_edep_%d",p+1),String.format("H_CND_t_edep_%d",p+1),100,0,250,100,T0,T1);
			H_CND_t_edep[p].setTitle(String.format("CND ver time vs E %d",p+1));
			H_CND_t_edep[p].setTitleX("Edep (MeV)");
			H_CND_t_edep[p].setTitleY("vt (ns)");
			H_CND_t[p] = new H1F(String.format("H_CND_t_%d",p+1),String.format("H_CND_t_%d",p+1),100,T0,T1);
			H_CND_t[p].setTitle(String.format("CND vertex time pad %d",p+1));
		}
		for(int iL=0;iL<3;iL++){
			H_CND_edep_z[iL] = new H2F("H_CND_edep_z","H_CND_edep_z",100,-50,70,100,0,250);
			H_CND_edep_z[iL].setTitle(String.format("CND Edep vs z",iL+1));
			H_CND_edep_z[iL].setTitleX("z (cm)");
			H_CND_edep_z[iL].setTitleY("E (MeV)");
			H_CND_edep_phi[iL] = new H2F("H_CND_edep_phi","H_CND_edep_phi",48,-180,180,100,0,250);
			H_CND_edep_phi[iL].setTitle(String.format("CND Edep vs #phi",iL+1));
			H_CND_edep_phi[iL].setTitleX("#phi (^o)");
			H_CND_edep_phi[iL].setTitleY("E (MeV)");
			H_CVT_CND_z[iL] = new H2F("H_CVT_CND_z","H_CVT_CND_z",100,-50,70,100,0,40);
			H_CVT_CND_z[iL].setTitle(String.format("CND z vs CVT z",iL+1));
			H_CVT_CND_z[iL].setTitleX("CND z (cm)");
			H_CVT_CND_z[iL].setTitleY("CVT z (cm)");
			H_CVT_CND_phi[iL] = new H2F("H_CVT_CND_phi","H_CVT_CND_phi",48,-180,180,100,-180,180);
			H_CVT_CND_phi[iL].setTitle(String.format("CND #phi vs CVT #phi",iL+1));
			H_CVT_CND_phi[iL].setTitleX("CND #phi (^o)");
			H_CVT_CND_phi[iL].setTitleY("CVT #phi (^o)");
			H_CND_vt_pad[iL] = new H2F(String.format("H_CND_phi_pad_"+(iL+1)),String.format("H_CND_phi_pad_"+(iL+1)),48,0.5,48.5,100,T0,T1);
			H_CND_vt_pad[iL].setTitle(String.format("CND vertex time vs pad",iL+1));
			H_CND_vt_pad[iL].setTitleX("CND counter #");
			H_CND_vt_pad[iL].setTitleY("CND vt (ns)");

			//H_neut_Edep_z, H_neut_z_phi, H_neut_bet_pad;
			//H_neut_Edep_z[iL] = new H2F("H_neut_Edep_z","H_neut_Edep_z",100,);
			//H_neut_Edep_z[iL].setTitle(String.format(""));
		}
	}
	public void FillCND(DataBank CNDbank, DataBank CVTbank){
		for(int iCND=0;iCND<CNDbank.rows();iCND++){
			int layer = CNDbank.getInt("layer",iCND);
			int trkID = CNDbank.getInt("trkID",iCND);
			if(layer>0 && layer<4 && trkID>-1){
				layer--;
				int sector = CNDbank.getInt("sector",iCND);
				int comp = CNDbank.getInt("component",iCND);
				float e  = CNDbank.getFloat("energy",iCND);
				float x  = CNDbank.getFloat("x",iCND);
				float y  = CNDbank.getFloat("y",iCND);
				float z  = CNDbank.getFloat("z",iCND);
				float tx = CNDbank.getFloat("tx",iCND);
				float ty = CNDbank.getFloat("ty",iCND);
				float tz = CNDbank.getFloat("tz",iCND);
				float time = CNDbank.getFloat("time",iCND);
				float path = CNDbank.getFloat("pathlength",iCND);
				float mom = CVTbank.getFloat("p",trkID);
				float beta = mom/(float)Math.sqrt(mom*mom+0.93827f*0.93827f);
				float vt = time - path/29.92f/beta - STT;
			
				float cndPhi = (float)Math.toDegrees(Math.atan2(y,x));
				float cvtPhi = (float)Math.toDegrees(Math.atan2(ty,tx));

				int pad = sector + 25*(comp-1);
				int PAD = pad + 50*layer;
				//vt -= Tmean[PAD-1];
				
				H_CND_occ.fill(pad);
				H_CND_time.fill(vt);
				if(PAD>0 && PAD<151){
					H_CND_t[PAD-1].fill(vt);
					H_CND_t_edep[PAD-1].fill(e,vt);
				}
				H_CND_phi_pad.fill(pad,cndPhi);
				H_CND_layer_pad.fill(layer+1 , pad);
				H_CND_vt_pad[layer].fill(pad,vt);
				H_CND_edep_z[layer].fill(z,e);
				H_CND_edep_phi[layer].fill(cndPhi,e);
				H_CVT_CND_z[layer].fill(z,tz);
				H_CVT_CND_phi[layer].fill(cndPhi,cvtPhi);
			}
		}
	}
	public void processEvent(DataEvent event) {
		DataBank eventBank = null;
		if(userTimeBased && event.hasBank("REC::Event"))eventBank = event.getBank("REC::Event");
		if(!userTimeBased && event.hasBank("RECHB::Event"))eventBank = event.getBank("RECHB::Event");
		
		if(eventBank!=null)STT = eventBank.getFloat("STTime",0);
		else return;
		if(event.hasBank("CND::hits") && event.hasBank("CVTRec::Tracks"))FillCND(event.getBank("CND::hits"),event.getBank("CVTRec::Tracks"));
	}
	public void plot() {
		EmbeddedCanvas can_cnd  = new EmbeddedCanvas();
                can_cnd.setSize(1500,3000);
		can_cnd.divide(3,6);
                can_cnd.setAxisTitleSize(18);
                can_cnd.setAxisFontSize(18);
                can_cnd.setTitleSize(18);
		can_cnd.cd(0);can_cnd.draw(H_CND_time);//can_cnd.draw(H_CND_occ);
		can_cnd.cd(1);can_cnd.draw(H_CND_phi_pad);
		can_cnd.cd(2);can_cnd.draw(H_CND_layer_pad);
		for(int iL=0;iL<3;iL++){
			can_cnd.cd(3+iL);can_cnd.draw(H_CND_vt_pad[iL]);
			can_cnd.cd(6+iL);can_cnd.draw(H_CND_edep_z[iL]);
			can_cnd.cd(9+iL);can_cnd.draw(H_CVT_CND_z[iL]);
			can_cnd.cd(12+iL);can_cnd.draw(H_CND_edep_phi[iL]);
			can_cnd.cd(15+iL);can_cnd.draw(H_CVT_CND_phi[iL]);
		}
                if(runNum>0){
                        if(!write_volatile)can_cnd.save(String.format("plots"+runNum+"/cnd.png"));
                        if(write_volatile)can_cnd.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/cnd.png"));
                        System.out.println(String.format("saved plots"+runNum+"/cnd.png"));
                }   
                else{
                        can_cnd.save(String.format("plots/cnd.png"));
                        System.out.println(String.format("saved plots/cnd.png"));
                }

		EmbeddedCanvas can_cnd2  = new EmbeddedCanvas();
                can_cnd2.setSize(10000,2400);
		can_cnd2.divide(25,6);
                can_cnd2.setAxisTitleSize(18);
                can_cnd2.setAxisFontSize(18);
                can_cnd2.setTitleSize(18);
		for(int p=0;p<150;p++){can_cnd2.cd(p);can_cnd2.draw(H_CND_t[p]);}
		if(runNum>0){
			if(!write_volatile)can_cnd2.save(String.format("plots"+runNum+"/cnd2.png"));
			if(write_volatile)can_cnd2.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/cnd2.png"));
			System.out.println(String.format("saved plots"+runNum+"/cnd2.png"));
		}else {
			can_cnd2.save("plots/cnd2.png");
			 System.out.println(String.format("saved plots/cnd2.png"));
		}

		EmbeddedCanvas can_cn3  = new EmbeddedCanvas();
                can_cn3.setSize(10000,2400);
		can_cn3.divide(25,6);
                can_cn3.setAxisTitleSize(18);
                can_cn3.setAxisFontSize(18);
                can_cn3.setTitleSize(18);
		for(int p=0;p<150;p++){can_cn3.cd(p);can_cn3.draw(H_CND_t_edep[p]);}
		if(runNum>0){
			if(!write_volatile)can_cn3.save(String.format("plots"+runNum+"/cnd3.png"));
			if(write_volatile)can_cn3.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/cnd3.png"));
			System.out.println(String.format("saved plots"+runNum+"/cnd3.png"));
		}
		else {
			can_cn3.save("plots/cnd3.png");
			System.out.println(String.format("saved plots/cnd3.png"));
		}


                for(int p=0;p<20;p++)System.out.print(String.format("%1.2ff , ",H_CND_t[p].getMean()));
                System.out.print("\n");
                for(int p=20;p<40;p++)System.out.print(String.format("%1.2ff , ",H_CND_t[p].getMean()));
                System.out.print("\n");
                for(int p=40;p<60;p++)System.out.print(String.format("%1.2ff , ",H_CND_t[p].getMean()));
                System.out.print("\n");
                for(int p=60;p<80;p++)System.out.print(String.format("%1.2ff , ",H_CND_t[p].getMean()));
                System.out.print("\n");
                for(int p=80;p<100;p++)System.out.print(String.format("%1.2ff , ",H_CND_t[p].getMean()));
                System.out.print("\n");
                for(int p=100;p<120;p++)System.out.print(String.format("%1.2ff , ",H_CND_t[p].getMean()));
                System.out.print("\n");
                for(int p=120;p<140;p++)System.out.print(String.format("%1.2ff , ",H_CND_t[p].getMean()));
                System.out.print("\n");
                for(int p=140;p<150;p++)System.out.print(String.format("%1.2ff , ",H_CND_t[p].getMean()));
                System.out.print("\n");
	}
        public void write(){
                TDirectory dirout = new TDirectory();
                dirout.mkdir("/cnd/");
                dirout.cd("/cnd/");
		for(int iL=0;iL<3;iL++)dirout.addDataSet(H_CND_vt_pad[iL],H_CND_edep_z[iL],H_CND_edep_phi[iL],H_CVT_CND_phi[iL]);
                for(int p=0;p<150;p++)dirout.addDataSet(H_CND_t[p]);
                if(write_volatile)if(runNum>0)dirout.writeFile("/volatile/clas12/rga/spring18/plots"+runNum+"/out_CND_"+runNum+".hipo");
                
		if(!write_volatile){
			if(runNum>0)dirout.writeFile("plots"+runNum+"/out_CND_"+runNum+".hipo");
			else dirout.writeFile("plots/out_CND.hipo");
		}
        }
////////////////////////////////////////////////
        public static void main(String[] args) {
                System.setProperty("java.awt.headless", "true");
                int count = 0;
		int runNum = 0;
		boolean useTB=true;
		boolean useVolatile = false;
		if(args.length>0)runNum=Integer.parseInt(args[0]);
		if(args.length>1)if(Integer.parseInt(args[1])==0)useTB=false;
		cnd ana = new cnd(runNum,useTB,useVolatile);
		HipoDataSource reader = new HipoDataSource();
		for(int ar=0;ar<args.length;ar++){
			reader.open(args[ar]);
			while(reader.hasEvent()) {
				DataEvent event = reader.getNextEvent();
				ana.processEvent(event);
				count++;
				if(count%10000 == 0) System.out.println(count/1000 + "k events (cnd on "+args[ar]+") ; "+(ar+1)+"/"+args.length);
			}
			reader.close();
		}
		System.out.println("Total : " + count + " events");
		ana.plot();
		ana.write();
        }
}
