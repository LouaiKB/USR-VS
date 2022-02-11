#### Overview

[USR-VS](http://usr.marseille.inserm.fr) is a web server application that implements the Ultrafast Shape Recognition algorithm, and uses two validated ligand-based 3D methods for large-scale prospective VS (USR and USRCAT). It's based on [istar's](https://jstar.cloud/) architecture which is for general-purpose web applications. This architecture is simplified with three main components: Web server, database server and daemon all deployed on a single virtual machine. In consequence, the network file system is replaced by an ordinary file system.

![Drag Racing](architecture-USR.png)

The daemon and core algorithms of this web server (USR and USRCAT) are written in C++14 for performance guarantee and memory management. An additional reason for using C++ over the other languages is the fact that [RDKit](https://www.rdkit.org/)  package provides a C++ API.
The web pages are built using Twitter’s [Bootstrap](https://getbootstrap.com/) template with HTML5, CSS3 & JavaScript. The web server is powered by [NodeJS](https://nodejs.org/en/), as it is fast and efficient and executes JavaScript.
The database server is powered by [MongoDB](https://www.mongodb.com/). It stores user-submitted jobs only, as ZINC molecules and conformers and their properties are stored in separate files on the file system.

All these components support Linux and Windows. Though this installation guide is for Linux, it is also possible to deploy USR-VS on a Windows server.


#### Installation

There is only one release of USR-VS and it was developed over six years. Therefore, the requirements of this web application are a bit old. We will also make some changes in some files to enable the safe installation.


#### Installing nodejs 0.10.40 for runing the web server

<br>

The web server of USR-VS was developed with the nodejs 0.10.40. Note that the prerequisites are for the target server of deployment. For instance, if USR-VS is deploy to a VM, these prerequisites are for the VM only, not the physical machine that hosts the VM. 
<br>

```{sh}
cd 
wget https://nodejs.org/dist/v0.10.40/node-v0.10.40-linux-x64.tar.gz
tar zxf node-v0.10.40-linux-x64.tar.gz && cd node-v0.10.40-linux-x64.tar.gz
sudo cp -rp node-v0.10.40-linux-x64 /usr/local/
sudo ln -s /usr/local/node-v0.10.40-linux-x64 /usr/local/node

```
<br>
Now we will add the node to the environment variables. For that open the `~/.bashrc`:
<br>
```
nano ~/.bashrc
```

Then add:


```
export PATH="/usr/local/node/bin/:$PATH"
```





USR
===

This webserver is based on [istar] and implements the ultrafast shape recognition algorithms.


Architecture
------------

![istar architecture](https://github.com/HongjianLi/istar/raw/master/public/architecture.png)


Components
----------

### Web client

* [jQuery]
* [Twitter Bootstrap]
* [three.js]
* [jquery-dateFormat]

### Web server

* [node.js]
* [mongodb]
* [express]

### Database

* [MongoDB]

### Daemon

* [MongoDB C++ Driver]


Supported browsers
------------------

* Mozilla Firefox 41


Licenses
--------

* Source code is licensed under [Apache License 2.0].
* Documentation is licensed under [CC BY 3.0].


Author
------

[Jacky Lee]


Logo
----

![USR logo](https://github.com/HongjianLi/usr/raw/master/public/logo.png)



[istar]: http://istar.cse.cuhk.edu.hk
[Twitter Bootstrap]: https://github.com/twitter/bootstrap
[jQuery]: https://github.com/jquery/jquery
[three.js]: https://github.com/mrdoob/three.js
[jquery-dateFormat]: https://github.com/phstc/jquery-dateFormat
[node.js]: https://github.com/joyent/node
[mongodb]: https://github.com/mongodb/node-mongodb-native
[express]: https://github.com/visionmedia/express
[MongoDB]: https://github.com/mongodb/mongo
[MongoDB C++ Driver]: https://github.com/mongodb/mongo-cxx-driver
[Apache License 2.0]: http://www.apache.org/licenses/LICENSE-2.0
[CC BY 3.0]: http://creativecommons.org/licenses/by/3.0
[Jacky Lee]: http://www.cse.cuhk.edu.hk/~hjli
