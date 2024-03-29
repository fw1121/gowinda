package gowinda.io;

import java.util.*;
import gowinda.analysis.GOEntry;
import java.io.*;
import gowinda.analysis.GOCategoryContainer;

/*
 * NOTE: Geneids are converted to lowercase
 */
public class GOEntryBulkReader {
	

	private String inputFile="";
	private java.util.logging.Logger logger;
	private BufferedReader bf;
	
	// Internal class containing the raw file entries
	private class GORaw
	{
		public String goname;
		public String description;
		public HashSet<String> geneids;
		public GORaw(String goname, String description, String[] geneids)
		{
			if(geneids.length <1 )throw new IllegalArgumentException("Your gene set file may be corrupt. At least one entry is required for every gene");
			this.goname=goname;
			this.description=description;
			
			this.geneids=new HashSet<String>();
			for(String s: geneids)
			{
				this.geneids.add(s.toLowerCase());
			}
		}

	}
	
	public GOEntryBulkReader(String inputFile, java.util.logging.Logger logger)
	{

		this.inputFile=inputFile;
		this.logger=logger;
		try
		{
			this.bf=new BufferedReader(new FileReader(inputFile));
		}
		catch(FileNotFoundException e)
		{
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	public GOEntryBulkReader(BufferedReader bf)
	{
		this.bf=bf;
		this.logger=java.util.logging.Logger.getLogger("Gowinda Logger");
		logger.setUseParentHandlers(false);
		
	}
	
	public GOCategoryContainer readGOEntries()
	{
		this.logger.info("Starting to read gene set file: "+this.inputFile);
		ArrayList<GORaw> ar= readGOraw();
		this.logger.info("Finished - read "+ar.size() + " gene sets");
		
		// Converting from GO: Geneid_list to Geneid vs. GOTerm_list 
		HashMap<String,HashSet<GOEntry>> goe=new  HashMap<String,HashSet<GOEntry>>();
		for(GORaw gr:ar)
		{
			GOEntry current=new GOEntry(gr.goname,gr.description);
			for(String geneid:gr.geneids)
			{
				if(!goe.containsKey(geneid))goe.put(geneid, new HashSet<GOEntry>());
				goe.get(geneid).add(current);
			}
		}
		
		HashMap<String,ArrayList<GOEntry>> toret=new HashMap<String,ArrayList<GOEntry>>();
		for(Map.Entry<String,HashSet<GOEntry>> en: goe.entrySet())
		{
			toret.put(en.getKey(), new ArrayList<GOEntry>(en.getValue()));
		}
		return new GOCategoryContainer(toret);
	}
	
	
	// Read the content of the GO file 
	private ArrayList<GORaw> readGOraw()
	{
		String line;
		ArrayList<GORaw> ar=new ArrayList<GORaw>();
		try
		{
			while((line=bf.readLine())!=null)
			{
				if(line.startsWith("#"))continue;
				String[] a=line.split("\t");
				String[] b={a[2]};
				if(a[2].contains(" ")) 
				{
					b=a[2].split(" ");
				}
				ar.add(new GORaw(a[0],a[1],b));
			}
			bf.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
			System.exit(0);
		}
		return ar;
	}
}