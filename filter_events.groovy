
import org.jlab.io.hipo.*;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

String inputfile = args[0];
int nevent = Integer.parseInt(args[1]);

HipoDataSync writer = new HipoDataSync();
writer.open("skim4.hipo");

int ntot=0;

HipoDataSource reader = new HipoDataSource();
reader.open(args[0]);

while(reader.hasEvent())
{

    if(ntot%10000 == 0) System.out.println("Analyzed " + ntot + " events");
    DataEvent event = reader.getNextEvent();
    // getting event number
    if (event.hasBank("RUN::config")) {
        DataBank bank = event.getBank("RUN::config");
       	ntot++;

        if(bank.getInt("event",0)==nevent) {
	   System.out.println("Found event " + nevent);
//	   event.show();
//           writer.writeEvent(event);
//           writer.writeEvent(event);
           writer.writeEvent(event);
           break;
        }
   }
}

writer.close();
