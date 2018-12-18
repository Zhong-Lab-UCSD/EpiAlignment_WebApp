const express = require('express')
const multer = require('multer')
const fs = require('fs');
const upload = multer({ dest: 'uploads/' })
const spawn = require('child_process').spawn
// express.js
const app = express()
// python
const pythonScript = 'server_agent.py'

var RunID_dict = {}
const RUNNING_CODE = -1

function makeid() {
  let text = ''
  let possible = 'abcdefghijklmnopqrstuvwxyz0123456789'

  for (var i = 0; i < 20; i++)
    text += possible.charAt(Math.floor(Math.random() * possible.length))

  return text
}

app.get('/', function (req, res) {
  res.sendFile( __dirname + '/' + 'index.html' )
})


const cpUpload = upload.fields([{ name: 'speciesPeak[]', maxCount: 2 }, { name: 'speciesInput[]', maxCount: 2 }])
app.post('/form_upload', cpUpload, function (req, res) {
  // req.files is an object (String -> Array) where fieldname is the key, and the value is array of files
  //
  // e.g.
  //  req.files['avatar'][0] -> File
  //  req.files['gallery'] -> Array
  //
  // req.body will contain the text fields, if there were any

  // Make new folder for results
  // flag = 1: folder exists
  var flag = 1
  do {
    try {
      var runid = makeid()
      tmp_stat = fs.statSync("tmp_" + runid)
    }
    catch(err){
      // if the folder does not exist
      console.log("Just created a new directory: " + runid)
      fs.mkdirSync("tmp_" + runid)
      RunID_dict[runid] = RUNNING_CODE
      flag = 0
    }
  }
  while (flag === 1)
 
  // construct an object for python input
  var pyMessenger = {'body': req.body, 'files': req.files, 'runid': "tmp_" + runid}
  var stdoutData = ''
  var stderrData = ''

  // execute python code on the server.
  let scriptExecution = spawn('python', [pythonScript]);

  // Handle normal output
  scriptExecution.stdout.on('data', (data) => {
      stdoutData += data
      console.log(data + '')
  })

  // Handle error output
  scriptExecution.stderr.on('data', (data) => {
    stderrData += data
    console.log(data + '')
  })

  scriptExecution.on('exit', (code) => {
      console.log("Process quit with code : " + code)
      console.log(stdoutData)
      console.log(stderrData)
      RunID_dict[runid] = code
  });
  // console.log(JSON.stringify(pyMessenger))
  // python input
  scriptExecution.stdin.write(JSON.stringify(pyMessenger));
  // tell the node that sending inputs to python is done.
  scriptExecution.stdin.end();

  res.json(req.files)

})

app.get('/download_results/:runid', function(req, res){
  var runid = req.params.runid
  console.log(runid)
  console.log(RunID_dict)
  if (!RunID_dict.hasOwnProperty(runid)) {
    // The query id does not exist.
    console.log(runid)
    res.send("This query ID does not exist!")
  } else if (RunID_dict[runid] === RUNNING_CODE) {
    res.send("EpiAlignment is still running.")
  } else if (RunID_dict[runid] === 0) {
    // python exited successfully.
    var file = __dirname + '/tmp_' + runid + '/' + 'epialign_result_' + runid
    console.log(file)
    res.sendFile(file);
  } else {
    // python exited with error.
    res.send("EpiAlignment exited with an error.")
  }
});

const server = app.listen(3000, "sysbio.ucsd.edu", function () {
 
  var host = server.address().address
  var port = server.address().port
 
  console.log('Listening http://%s:%s', host, port)
 
})