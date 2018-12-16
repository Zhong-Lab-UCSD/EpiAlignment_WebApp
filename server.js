const express = require('express')
const multer = require('multer')
const fs = require('fs');
const upload = multer({ dest: 'uploads/' })
const spawn = require('child_process').spawn
// express.js
const app = express()
// python
const pythonScript = 'server_agent.py'

var stdoutData = ''
var stderrData = ''

function makeid() {
  let text = ''
  let possible = 'abcdefghijklmnopqrstuvwxyz0123456789'

  for (var i = 0; i < 20; i++)
    text += possible.charAt(Math.floor(Math.random() * possible.length))

  return text
}

app.get('/', function (req, res) {
  console.log(__dirname + '/' + 'index.html')
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
      flag = 0
    }
  }
  while (flag === 1)
 
  // construct an object for python input
  var pyMessenger = {'body': req.body, 'files': req.files, 'runid': "tmp_" + runid}

  // execute python code on the server.
  let scriptExecution = spawn('python', [pythonScript]);

  // Handle normal output
  scriptExecution.stdout.on('data', (data) => {
      stdoutData += data
  })

  // Handle error output
  scriptExecution.stderr.on('data', (data) => {
    stderrData += data
  })

  scriptExecution.on('exit', (code) => {
      console.log("Process quit with code : " + code)
      console.log(stdoutData)
      console.log(stderrData)
  });
  console.log(JSON.stringify(pyMessenger))
  // python input
  scriptExecution.stdin.write(JSON.stringify(pyMessenger));
  // tell the node that sending inputs to python is done.
  scriptExecution.stdin.end();

  res.json(req.files)
  res.end()

})

const server = app.listen(3000, function () {
 
  var host = server.address().address
  var port = server.address().port
 
  console.log('Listening http://%s:%s', host, port)
 
})