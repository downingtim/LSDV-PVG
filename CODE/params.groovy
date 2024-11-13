// Initialize an empty map to hold the parameters
def params = [:]
args.each { arg ->
    def (key, value) = arg.split('=')
    params[key] = value
}

for(Lib in Keys)
{
    Options = YAML.get(Lib)
    
    def getParam = { param, optParam -> param ?: (optParam ?: "") }
    
    if (Executor_Parameters) {
        Q_List.push(getParam(Executor_Parameters.queue, Options?.process?.queue))
        Memory.push(getParam(Executor_Parameters.memory, Options?.process?.memory))
        Cpus_List.push(getParam(Executor_Parameters.cpus, Options?.process?.cpus))
        ClusterOptions_List.push(getParam(Executor_Parameters.clusterOptions, Options?.process?.clusterOptions))
        Executor = Executor_Parameters.executor ?: ""
    } else {
        Q_List.push(Options?.process?.queue ?: "")
        Memory.push(Options?.process?.memory ?: "")
        Cpus_List.push(Options?.process?.cpus ?: "")
        ClusterOptions_List.push(Options?.process?.clusterOptions ?: "")
    }
    
    if (!Executor_Parameters && !Options?.parameters) {
        Q_List.push("");Memory.push("");Cpus_List.push("");ClusterOptions_List.push("")
    }
}


/*static def Make_Dir(String S) 
{
  if(!file(S).mkdirs())
  {
    println ("Make_Dir : Cannot create directory " + S);
    exit(-1);
  }
}*/


/*if(params.find { it.key == 'config'} == null )
{ println "Please specify a configuration file with the option --config \n";
  exit(-1); 
}
if(params.find { it.key == 'trimmer'} == null )
{
  params.trimmer = 1; 
}
if(params.find { it.key == 'assembler'} == null )
{
  params.assembler = 1; 
}
if(params.find { it.key == 'finisher'} == null )
{
  params.finisher = 1; 
}
if(params.find { it.key == 'variantcaller'} == null )
{
  params.variantcaller = 1; 
}
if(params.find { it.key == 'stat'} == null )
{
  params.stat = 1; 
}
if(params.find { it.key == 'identity'} == null )
{
  params.identity = 80; 
}
*/


//Load YAML configuration..  
@Grab('org.yaml:snakeyaml:1.17')

import org.yaml.snakeyaml.Yaml

Yaml parser = new Yaml()
YAML = parser.load((params.config as File).text);
Keys = YAML.keySet();

//Set Configuration ..
Executor_Parameters = null;
if(Keys.find{ it == 'PROCESS' } != null )//Process specific parameters..
{
  Executor_Parameters = YAML.get('PROCESS');
  Keys.removeAll{ it == 'PROCESS'};
}

Q_List = [];
Memory = [];
Cpus_List = [];
ClusterOptions_List = [];
Executor = "";

/*
for(Lib in Keys)
{
  Options = YAML.get(Lib);
  if(Executor_Parameters != null)
  {
    if(Options.parameters != null)
    {
      Q_List.push(( Executor_Parameters.queue || Options.process.queue ) ? ((Options.process.queue) ? Options.process.queue : Executor_Parameters.queue) : "" ); 
      Memory.push(( Executor_Parameters.memory  || Options.process.memory) ?  ((Options.process.memory) ? Options.process.memory : Executor_Parameters.memory) : "") ; 
      Cpus_List.push(( Executor_Parameters.cpus || Options.process.cpus ) ?  ((Options.process.cpus) ? Options.process.cpus : Executor_Parameters.cpus) : "") ; 
      ClusterOptions_List.push(( Executor_Parameters.clusterOptions || Options.process.clusterOptions ) ? ((Options.process.clusterOptions) ? Options.process.clusterOptions : Executor_Parameters.clusterOptions) : ""); 
      Executor =( Executor_Parameters.executor ) ? Executor_Parameters.executor : "";
    }
    else
    {
      Q_List.push(( Executor_Parameters.queue ) ? Executor_Parameters.queue : "" ) ; 
      Memory.push(( Executor_Parameters.memory ) ?  Executor_Parameters.memory : "") ; 
      Cpus_List.push(( Executor_Parameters.cpus ) ?  Executor_Parameters.cpus : "") ; 
      ClusterOptions_List.push(( Executor_Parameters.clusterOptions ) ? Executor_Parameters.clusterOptions : ""); 
      Executor =( Executor_Parameters.executor ) ? Executor_Parameters.executor : "";
    }
  }
  else
  {
    if(Options.parameters != null)
    {
      Q_List.push((Options.process.queue) ? Options.process.queue : "" ); 
      Memory.push((Options.process.memory) ? Options.process.memory : "") ; 
      Cpus_List.push((Options.process.cpus) ? Options.process.cpus : "") ; 
      ClusterOptions_List.push((Options.process.clusterOptions) ? Options.process.clusterOptions : ""); 
    }
    else
    {
      Q_List.push( "" ); Memory.push( "" ); Cpus_List.push( "" ); ClusterOptions_List.push( "" );  
    }
  }
}

//Create channels ..
Lib_List = Channel.from(Keys);*/
