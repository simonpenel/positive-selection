var express = require('express');
var router = express.Router();
const browser = require('browser-detect');
var bodyParser = require('body-parser');
router.use(bodyParser.json()); // support json encoded bodies
router.use(bodyParser.urlencoded({ extended: true })); // support encoded bodies
if (typeof localStorage === "undefined" || localStorage === null) {
  var LocalStorage = require('node-localstorage').LocalStorage;
  localStorage = new LocalStorage('./scratch');
}
const util = require('util')

// =====================
// Routages GET
// =====================

// Home
// ----
router.get('/', function(req, res, next) {
  res.render('index', { title: 'HOGENOM'});
});

// GET display
// -----------
router.get('/display', function(req, res, next) {
  const fs = require('fs');
  fs.readFile('input_tree.xml', 'utf8' , (err, data) => {
    if (err) {
      res.render('error.ejs', {message:"Erreur de lecture",error:err});
    }
    var xml_digester = require("xml-digester");
    var handler = new xml_digester.OrderedElementsHandler("eventType");
    var options = {
      "handler": [{
        "path": "eventsRec/*",
        "handler": handler
      }]
    };
    var digester = xml_digester.XmlDigester(options);
    digester.digest(data, function(err, results) {
      if (err) {
        console.log(err);
        return
      }
      var JSONtree = JSON.stringify(results);
      var JSONpattern = JSON.stringify("0:DELAS_1.PE4738");
      res.render('displaytree.ejs', {arbre:JSONtree,pattern:JSONpattern});
    });
  });
});

// GET taxodico
// ------------
router.get('/taxodico', function(req, res, next) {
  res.render('gettaxojson.ejs', {species:JSON.stringify([]),colour:JSON.stringify([])});
});

module.exports = router;
