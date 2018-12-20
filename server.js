const express = require('express')
const multer = require('multer')
const fs = require('fs');
const nodemailer = require('nodemailer')
const upload = multer({ dest: 'uploads/' })
const spawn = require('child_process').spawn
// express.js
const app = express()
// python
const pythonScript = 'server_agent.py'
// util
const util = require('util')
const readFilePromise = util.promisify(fs.readFile)

var RunID_dict = {}
const RUNNING_CODE = -1

function makeid() {
  let text = ''
  let possible = 'abcdefghijklmnopqrstuvwxyz0123456789'

  for (var i = 0; i < 20; i++)
    text += possible.charAt(Math.floor(Math.random() * possible.length))

  return text
}

const cpUpload = upload.fields([{ name: 'speciesPeak1[]', maxCount: 3 }, { name: 'speciesPeak2[]', maxCount: 3 }, { name: 'speciesInput1', maxCount: 1 }, { name: 'speciesInput2', maxCount: 1 }])

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
    catch (err) {
      // if the folder does not exist
      console.log("Just created a new directory: " + runid)
      fs.mkdirSync("tmp_" + runid)
      RunID_dict[runid] = RUNNING_CODE
      flag = 0
    }
  }
  while (flag === 1)

  console.log(RunID_dict)

  // construct an object for python input
  var pyMessenger = { 'body': req.body, 'files': req.files, 'runid': "tmp_" + runid }
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
    let message = null
    if (req.body.emailNote && req.body.mail) {
      message = {
        from: 'no-reply@mail.beta.givengine.org',
        to: req.body.mail
      }
    }
    let header = 'Dear user,\n\n'
    let footer = '\n\nThe EpiAlignment Team\n--\nZhong Lab\n' +
      'Department of Bioengineering, Mail Code 0412\n' +
      'University of California, San Diego\n' +
      '9500 Gilman Dr.\nLa Jolla, CA 92122-0412\nUnited States'

    if (message) {
      console.log("Process quit with code : " + code)
      RunID_dict[runid] = code
      // needs to write an email
      let transporter = nodemailer.createTransport({
        sendmail: true,
        newline: 'unix',
        path: '/usr/sbin/sendmail'
      })
      if (code) {
        // some kind of error has happened
        message.subject = 'Error when computing your EpiAlignment data'
        message.text = header +
          'Unfortunately, an error occured when processing ' +
          'your data for EpiAlignment (Error Code: ' + code +
          '). We are unable to complete your request. \n\n' +
          'In some cases this may be due to an erroneous input format, ' +
          'in which case you may try again by providing the correctly ' +
          'formatted input at https://beta.epialign.ucsd.edu/.\n\n' +
          'If the error keep happening, please let us know by emailing ' +
          'Jia Lu<jil430@eng.ucsd.edu> or Xiaoyi Cao <x9cao@eng.ucsd.edu>. ' +
          'Thank you for trying EpiAlignment!' + footer
      } else {
        // Normal email
        message.subject = 'Your EpiAlignment results are ready'
        message.text = header +
          'Your EpiAlignment results are now ready and can be downloaded ' +
          'by following this link:\n\n' +
          'https://beta.epialign.ucsd.edu/result_page/' + runid + '\n\n' +
          'If the link above does not work, please copy the entire link ' +
          'and paste it into the address bar of your web browser.\n\n' +
          'If you have any questions, please let us know by emailing ' +
          'Jia Lu<jil430@eng.ucsd.edu> or Xiaoyi Cao <x9cao@eng.ucsd.edu>.' +
          ' \n\nThank you for using EpiAlignment!' + footer
      }
      transporter.sendMail(message)
    }
  });
  // console.log(JSON.stringify(pyMessenger))
  // python input
  scriptExecution.stdin.write(JSON.stringify(pyMessenger));
  // tell the node that sending inputs to python is done.
  scriptExecution.stdin.end();

  console.log(RunID_dict)
  res.json({ runid: runid })

})

app.get('/results/:runid', function (req, res) {
  var runid = req.params.runid
  console.log(runid)
  console.log(RunID_dict)
  if (!RunID_dict.hasOwnProperty(runid)) {
    // The query id does not exist.
    console.log(runid)
    res.status(401)
    res.send("This query ID does not exist!")
  } else if (RunID_dict[runid] === RUNNING_CODE) {
    res.json({ status: RUNNING_CODE })
  } else if (RunID_dict[runid] === 0) {
    // python exited successfully.
    // read JSON file and return it.
    readFilePromise('tmp_' + runid + '/' + runid + '.json', 'utf8')
      .then(result => {
        let resObj = JSON.parse(result)
        res.json({
          status: 0,
          data: resObj
        })
      })
      .catch(err => {
        res.status(401)
        res.send(err.message)
      })
  } else {
    // python exited with error.
    res.json({
      // TODO: find a way to insert error code and msg here
      status: 1,
      message: "EpiAlignment exited with an error."
    })
  }
});

const server = app.listen(3000, function () {

  var host = server.address().address
  var port = server.address().port

  console.log('Listening http://%s:%s', host, port)

})